# Load the libraries and set working directory
source("analysis/setup.R")

## load region boundary
cfr=readOGR("data/src/CFR","CFR")
ereg=extent(cfr)

cf=stack(c("data/MCD09/MCD09_mean_01.tif","data/MCD09/MCD09_mean_07.tif","data/MCD09_deriv/intra.tif"))
names(cf)=c("CLDJAN","CLDJUL","CLDSEAS")
NAvalue(cf)=0
cf=crop(cf,cfr)
cf=mask(cf,cfr)
gain(cf)=.01

sp="PRCYNA"
sp2="Protea_cynaroides"

## create output folder
outputdir=paste0("output/sdm/",sp2,"/")
if(!file.exists(outputdir)) dir.create(outputdir,recursive=T)

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
env=mask(env,cfr)

cmeans=cellStats(env,"mean")
csd=cellStats(env,"sd")
senv=scale(env)


## rasterize protea data to grid and return a raster of presences and trials
fprot=function(prot,rast=cf,sp="PRACAU"){
  tpres=unique(prot@data[prot$pro==sp,c("obs","londd","latdd")])
  coordinates(tpres)=c("londd","latdd")
  ttrials=unique(prot@data[,c("obs","londd","latdd")])
  coordinates(ttrials)=c("londd","latdd")
  presences=rasterize(tpres,cf,fun="count",field="obs",background=0)
  trials=rasterize(ttrials,cf,fun="count",field="obs",background=0)
  td=mask(stack(presences,trials),cfr)
  names(td)=c("presences","trials")
  return(td)
}

rprot=fprot(prot,sp=sp)

## build big raster stack
fit=rprot[["trials"]]>0
rdata=stack(senv,rprot,fit)
writeRaster(rdata,file=paste0(outputdir,"/",sp2,"_modelinput.nc"),overwrite=T)#,coordinates(senv),cell=cellFromXY(senv,xy=coordinates(senv)))

system(paste0("ncdump -h ",paste0(outputdir,"/",sp2,"_modelinput.nc")))

data=cbind.data.frame(values(stack(senv,rprot)),coordinates(senv),cell=cellFromXY(senv,xy=coordinates(senv)))
data=as.data.frame(data[!is.na(data[,"presences"]),])
data$trials=as.integer(data$trials)
data$presences=as.integer(data$presences)
data=na.omit(data)
data$fit=ifelse(data$trials>0,1,0)



## create 'fitting' dataset where there are observations
fdata=data[data$fit==1,]


### Save model input data
save(env,senv,data,fdata,sp,sp2,file=paste0(outputdir,"modelinput.Rdata"))

##############################################################
load(paste0(outputdir,"modelinput.Rdata"))

nchains=3

mods=data.frame(
  model=c("m1","m2"),
  formula=c("~PPTJAN+PPTJUL+PPTSEAS+MAT",
            "~CLDJAN+CLDJUL+CLDSEAS+MAT"),
  name=c("Precipitation",
         "Cloud"))


registerDoMC(ncores)


foreach(m=1:nrow(mods)) %dopar% { 
      fsdm(fdata=fdata,data=data,model=mods$formula[m],modelname=mods$name[m],nchains=nchains,
           species=sp2,outputdir=outputdir,verbose=T,autocor=T)    
    }
    

## regional
p1=ggplot(pred) + geom_tile(aes(x=x,y=y,fill = pred)) +
  facet_wrap(~modelname,ncol=1) +
  scale_fill_gradientn(colours=c('white','blue','red'),
                       name="P(Presence)") +
  coord_equal()+
  xlim(c(20,25))+ylim(-34.2,-33.5)+
#  geom_point(data=prot[prot$pro!=sp,]@data,aes(x=londd,y=latdd),col="grey",alpha=.2,pch=1)+
#  geom_point(data=prot[prot$pro==sp,]@data,aes(x=londd,y=latdd),pch="+",cex=1)+
  theme(panel.background = element_rect(fill='transparent'),legend.key.width = unit(1.5, "cm"))+
  xlab(label="Longitude")+ylab("Latitude")