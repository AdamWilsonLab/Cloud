# Load the libraries and set working directory
source("analysis/setup.R")

sp="PRCYNA"
sp2="Protea_cynaroides"

## create output folder
outputdir=paste0("output/sdm/",sp2,"/")
if(!file.exists(outputdir)) dir.create(outputdir,recursive=T)
fmodelinput=paste0(outputdir,"/",sp2,"_modelinput.nc")

## load region boundary
cfr=readOGR("data/src/CFR","CFR")
ereg=extent(cfr)

cf=stack(c("data/MCD09/MCD09_mean_01.tif","data/MCD09/MCD09_mean_07.tif","data/MCD09_deriv/intra.tif"))
names(cf)=c("CLDJAN","CLDJUL","CLDSEAS")
NAvalue(cf)=0
cf=crop(cf,cfr)
gain(cf)=.01


## Load protea data
prot=readOGR("data/src/proteaatlas/","prot")
prot=prot[prot$pla%in%c("NA","A","E"),]  #get rid of planted

prot=prot[-grep(".",prot$pro,fixed=T),]  # get rid of hybrids with .'s
prot$pro=as.factor(as.character(prot$pro)) #reset the factors


### Load environmental data
tvars=c(
  paste0("/mnt/data/jetzlab/Data/environ/global/worldclim/",
         c("prec_1","prec_7","bio_1","bio_12","bio_15","tmean_1","tmean_7","alt"),".bil"))
tenv=stack(crop(stack(tvars),ereg),cf)

## Covariate correlation
tcor=layerStats(tenv, "pearson",asSample=F, na.rm=T)
write.csv(tcor,file=paste0(outputdir,sp2,"_covariatecorrelation.csv"),row.names=F)

### Select environmental data
##  Create single stack
env=stack(tenv[[c("prec_1","prec_7","bio_15","bio_1","alt")]],cf[[c("CLDJAN","CLDJUL","CLDSEAS")]])
names(env)=c("PPTJAN","PPTJUL","PPTSEAS","MAT","ALT","CLDJAN","CLDJUL","CLDSEAS")

env=crop(env,ereg)
#env=mask(env,cfr)


fprot=function(prot,rast=cf,sp="PRACAU"){
  obs=unique(prot@data[,c("obs","londd","latdd")])
  obs$trial=1
  obs$presence=ifelse(obs$obs%in%prot$obs[prot$pro==sp],1,0)
  coordinates(obs)=c("londd","latdd")
  obs@data[,c("lon","lat")]=coordinates(obs)
  return(obs)
}

pprot=fprot(prot,sp=sp)



## define metadata that will be embeded in output object
meta=list(institution="Map of Life, Yale University, New Haven, CT",
          source="Species Distributions",
          comment="Adam M. Wilson (adam.wilson@yale.edu)",
          species=sp2)

hSDM.ncWriteInput(env,points=pprot,ncfile=fmodelinput,overwrite=T,meta=meta)


##############################################################
## import data
data=hSDM.ncReadInput(fmodelinput)
## select data for fitting
data$fit=ifelse(data$trials>0,1,0)
## create 'fitting' dataset where there are observations
fdata=data[data$fit>0,]


nchains=3
registerDoMC(ncores)

mods=data.frame(
  model=c("m1","m2"),
  formula=c("~PPTJAN+PPTJUL+PPTSEAS+MAT",
            "~CLDJAN+CLDJUL+CLDSEAS+MAT"),
  name=c("Precipitation",
         "Cloud"))


registerDoMC(ncores)

burnin=10000
mcmc=10000
thin=10

foreach(m=1:nrow(mods)) %dopar% { 
  tm=P.hSDM.ZIB(
    suitability=as.character(mods$formula[m]),
    presences=fdata$presences,
    observability=~1,
    mugamma=0, Vgamma=1.0E6,
    gamma.start=0,
    trials=fdata$trials,
    data=fdata,
    burnin=burnin, mcmc=mcmc, thin=thin,
    beta.start=0,
    suitability.pred=data,
    mubeta=0, Vbeta=1.0E6,
    save.p=0,
    verbose=1,
    seed=round(runif(1,0,1e6)))

  meta=list(modelname=as.character(mods$name[m]),species=sp2)
  outfile=paste0(outputdir,"/",sp2,"_modeloutput_",mods$name[m],".nc")
 source("/media/data/repos/hsdm-code/R/hSDM.ncWriteOutput.R")     
  hSDM.ncWriteOutput(results=tm,file=outfile,meta=meta,verbose=T,autocor=T)    
    }


## regional
#  xlim(c(20,25))+ylim(-34.2,-33.5)+
