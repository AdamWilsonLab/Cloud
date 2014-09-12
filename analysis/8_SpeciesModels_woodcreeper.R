# Load the libraries and set working directory
source("analysis/setup.R")

registerDoMC(12)
rasterOptions(datatype="FLT4S")
  

sp=c("Lepidocolaptes","lacrymiger")
fam="Furnariidae (Ovenbirds and Woodcreepers)"

## load region boundary
download.file(paste0("http://mol.cartodb.com/api/v2/sql?",
                  "q=SELECT%20ST_TRANSFORM(the_geom_webmercator,4326)%20as%20the_geom,%20seasonality%20FROM%20",
                  "get_tile('jetz','range','",
                  sp[1],"%20",sp[2],
                  "','jetz_maps')&format=shp&filename=expert"),
              destfile="data/SDM/expert.zip")
unzip("data/SDM/expert.zip",exdir="data/SDM/")

reg=readOGR("data/SDM/","expert")

### Load eBird data for 
## ebird data downloaded and unzipped from http://ebirddata.ornith.cornell.edu/downloads/erd/ebird_all_species/erd_western_hemisphere_data_grouped_by_year_v5.0.tar.gz
ebirddir="/mnt/data/jetzlab/Data/specdist/global/ebird/"

## select species for presences and 'absences'
taxon=read.csv(paste0(ebirddir,"/doc/taxonomy.csv"))
taxon$TAXON_ORDER2=round(taxon$TAXON_ORDER)
taxon=taxon%.%filter(FAMILY_NAME==fam)

library(data.table)
library(AUC)
library(dismo)



fspdata=paste0("data/SDM/points_",paste(sp,collapse="_"),".csv")
if(!file.exists(fspdata)){
    f="/mnt/data2/projects/mol/points_Aug_2014/ebd_relMay-2014.txt"
    hdr=colnames( read.table(f,nrows=3,header=T,sep="\t"))
  ## define columns to keep  
    cols=c("GLOBAL UNIQUE IDENTIFIER","TAXONOMIC ORDER","CATEGORY",
         "COMMON NAME","SCIENTIFIC NAME",
         "OBSERVATION COUNT","COUNTRY","COUNTRY CODE",  
         "LATITUDE","LONGITUDE","OBSERVATION DATE",
         "OBSERVER ID","SAMPLING EVENT.IDENTIFIER","PROTOCOL TYPE",
         "DURATION MINUTES","EFFORT DISTANCE KM","EFFORT AREA HA",
         "NUMBER OBSERVERS","ALL SPECIES REPORTED","GROUP IDENTIFIER",
         "APPROVED","REVIEWED","REASON")

  # make a new copy with no quotes
  system(paste("sed -e 's|[\"'\']||g' ",f," > ebird.txt"))
    # 50 minutes
  d=fread("ebird.txt",header = T, sep = '\t',select=cols,na.strings=c("","?","\""),showProgress=T,verbose=T)
  setnames(d,colnames(d),gsub(" ","_",colnames(d)))
  setkey(d,"TAXONOMIC_ORDER")
    # subset by lat/lon - 2 minutes
  system.time(d2<<-d%.%mutate(TAXO=round(TAXONOMIC_ORDER))%.%
                filter(LATITUDE>=bbox(reg)["y","min"],
                       LATITUDE<=bbox(reg)["y","max"],
                       LONGITUDE>=bbox(reg)["x","min"],
                       LONGITUDE<=bbox(reg)["x","max"],
                       TAXO%in%taxon$TAXON_ORDER2))  
     write.csv(d2,file=fspdata,row.names=F)
  ## clean up
    rm(d,d2); gc()
    file.remove("ebird.txt")
}

####################
## read it back in


#########  Environmental Data
cf=stack(c("data/MCD09/MCD09_mean_01.tif","data/MCD09/MCD09_mean_07.tif","data/MCD09_deriv/meanannual.tif","data/MCD09_deriv/intra.tif"))
names(cf)=c("c_mm_01","c_mm_07","c_ma","c_intra")
NAvalue(cf)=0
cf=crop(cf,reg)
gain(cf)=.01

### Load environmental data
vars=c(
  paste0("/mnt/data/jetzlab/Data//environ/global/worldclim/",c("tmean_1","tmean_7","prec_1","prec_7","alt"),".bil"))
env=stack(vars)
names(env)=gsub("_","",names(env))
env=crop(env,reg)
##  Create single stack
env=stack(env,cf)


### clean up species data
spd=read.csv(fspdata)

sptax=round(taxon$TAXON_ORDER[taxon$GENUS_NAME==sp[1]&taxon$SPECIES_NAME==sp[2]])

## observations of all species
spd_ns=spd%.%filter(ALL_SPECIES_REPORTED==1)%.%
  group_by(LATITUDE,LONGITUDE,OBSERVATION_DATE)
spd_ns=unique(spd_ns[,c("LATITUDE","LONGITUDE","OBSERVATION_DATE")])%.%mutate(presence=0)

## presences
spd_s=spd%.%filter(round(TAXONOMIC_ORDER)==sptax)%.%select(LONGITUDE,LATITUDE,OBSERVATION_DATE)%.%mutate(presence=1)

spd_sp=as.data.frame(rbind(spd_ns,spd_s))
coordinates(spd_sp)=c("LONGITUDE","LATITUDE")
spd_sp@data[,c("lon","lat")]=coordinates(spd_sp)  
  
presences=rasterize(spd_sp[spd_sp$presence==1,],cf,fun="count",field="presence",background=0)
  trials=rasterize(spd_sp[spd_sp$presence==0,],cf,fun="count",field="presence",background=0)
  td=stack(presences,trials)
  names(td)=c("presences","trials")

senv=scale(env)

#splom(senv)

data=cbind.data.frame(values(stack(senv,td)),coordinates(senv),cell=cellFromXY(senv,xy=coordinates(senv)))
data=as.data.frame(data[!is.na(data[,"presences"]),])
data$trials=as.integer(data$trials)
data$presences=as.integer(data$presences)
data=na.omit(data)

## adjacency matrix
#neighbors.mat <- adjacent(senv, cells=1:ncell(senv), directions=8, pairs=TRUE, sorted=TRUE)
#neighbors.mat=neighbors.mat[neighbors.mat[,1]%in%data$cell,]
#n.neighbors <- as.data.frame(table(as.factor(neighbors.mat[,1])))[,2]
#adj <- neighbors.mat[,2]

## create 'fitting' dataset where there are observations
fdata=data[data$trials>0,]
## due to opportunistic observations, there are a few sites with more presences than trials, update those records so presences=trials
fdata$trials[fdata$presences>fdata$trials]=fdata$presences[fdata$presences>fdata$trials]

nchains=3

mods=data.frame(
  model=c("m1","m2","m3"),
  formula=c("~tmean1+tmean7+prec1+prec7+alt",
                "~tmean1+tmean7+c_mm_01+c_mm_07+alt",
                "~tmean1+tmean7+prec1+prec7+c_mm_01+c_mm_07+alt"),
  name=c("Temperature & Precipitation",
         "Temperature & Cloud",
         "Temperature, Precipitation, & Cloud"))

## file to save output
fres=paste0("data/SDM/",paste(sp,collapse="_"),"_modeloutput.Rdata")

if(!file.exists(fres)){
res=foreach(m=1:nrow(mods)) %dopar% { 
  tres1=foreach(ch=1:nchains) %dopar% {
    tres2=hSDM.ZIB(
      suitability=as.character(mods$formula[m]),
      presences=fdata$presences,
      observability=~1,
      mugamma=0, Vgamma=1.0E6,
      gamma.start=0,
      trials=fdata$trials,
      data=fdata,
      burnin=10000, mcmc=10000, thin=10,
      beta.start=0,
      suitability.pred=data,
      mubeta=0, Vbeta=1.0E6,
      save.p=0,
      verbose=1,
      seed=round(runif(1,0,1000)))
}

## combine the chains for each model
  coef=data.frame(model=mods$model[m],
                   model=mods$model[m],
                  modelname=mods$name[m],
                   param=colnames(tres1[[1]][[1]]),
                  summary(as.mcmc.list(lapply(tres1,FUN=function(x) x$mcmc)))$statistics,
                  summary(as.mcmc.list(lapply(tres1,FUN=function(x) x$mcmc)))$quantiles)
  pred=data.frame(model=mods$model[m],
                  model=mods$model[m],
                  modelname=mods$name[m],
                  x=data$x,y=data$y,
                  pred=rowMeans(do.call(cbind,lapply(tres1,FUN=function(x) x$prob.p.pred))))
  if(!is.null(tres1[[1]]$rho.pred)){
    rho=data.frame(model=mods$model[m],
                  suitability=mods$suitability[m],
                   model=mods$model[m],
                   modelname=mods$name[m],
                   coordinates(senv),
                   cell=1:ncell(senv),
                  rho=rowMeans(do.call(cbind,lapply(tres1,FUN=function(x) x$rho.pred))))
    rho=rho[rho$cell%in%data$cell,]
  }
  if(is.null(tres1[[1]]$rho.pred)){
    rho=NULL
  }
return(list(coef=coef,pred=pred,rho=rho))  
}


save(res,file=fres)

}

#################
load(fres)

coef=do.call(rbind,lapply(res,function(x) x$coef))
pred=do.call(rbind,lapply(res,function(x) x$pred))
rho=do.call(rbind,lapply(res,function(x) x$rho))

pred=pred[pred$model%in%mods$model[1:2],]
pred$cell=cellFromXY(senv,xy=pred[,c("x","y")])

predr=stack(lapply(unique(pred$model),function(m) rasterFromXYZ(xyz=pred[pred$model==m,c("x","y","pred")])))
names(predr)=unique(pred$model)
projection(predr)='+proj=longlat'
#writeRaster(predr*100,file=paste0("data/SDM/",paste(sp,collapse="_"),".tif"))


##### Autocorrelation
library(devtools)
install_github("adammwilson/rasterAutocorr")
library(rasterAutocorr)

ac1=acor_table(predr[[1]],verbose=T)
ac1$model=names(predr)[1]
ac2=acor_table(predr[[2]],verbose=T)
ac2$model=names(predr)[2]

ac=rbind.data.frame(ac1,ac2)


## Global Spatial Autocorrelation
fspac=paste0("output/sdm/globalautocorr_",paste(sp,collapse="_"),".csv")
if(!file.exists(fspac)){
spac=do.call(rbind.data.frame,
             lapply(1:nlayers(predr),function(i) {
                data.frame(species=paste(sp,collapse=" "),
                           model=names(predr)[i],
                           MoransI=Moran(predr[[i]],w=matrix(1,51,51)),
                           GearyC=Geary(predr[[i]],w=matrix(1,51,51)))}
))
write.csv(spac,file=fspac,row.names=F)
}
spac=read.csv(fspac)


## Range size
rs=cellStats(predr,"sum")
(rs[[1]]-rs[[2]])/rs[[1]]

splom(predr)

cor(fdata$c_mm_07,fdata$prec7)

# http://www.bayesian-inference.com/modelfit
# http://voteview.com/DIC-slides.pdf
dic=coef[grepl("Deviance",coef$param),]
dic$pv=(dic$SD^2)/2
dic$dic=dic$Mean+dic$pv

## AUC
aucdat=merge(fdata[,c("presences","trials","cell")],pred,by=c("cell"))
aucdat$model=factor(aucdat$model,levels=c("m1","m2"),ordered=T)

feval=function(x){
  e=evaluate(p=x$pred[x$presences>0],a=x$pred[x$presences==0])
  data.frame(nPresence=e@np,nTrials=e@na,auc=e@auc,cor=e@cor)
}

restable=do.call(rbind.data.frame,by(aucdat,INDICES=aucdat$modelname,FUN=feval))
restable$species=paste(sp,collapse=" ")
restable$modelname=unique(as.character(aucdat$modelname))
restable$model=dic$model[match(restable$modelname,dic$modelname)]
restable$Deviance=dic$Mean[match(restable$model,dic$modelname)]
restable$Pv=dic$pv[match(restable$model,dic$modelname)]
restable$DIC=dic$dic[match(restable$model,dic$modelname)]
##
restable[,c("MoransI","GearyC")]=spac[match(restable$model,spac$model),c("MoransI","GearyC")]

restable=restable[,c("species","model","modelname","nPresence","nTrials","Deviance","Pv","DIC","auc","cor","MoransI","GearyC")]

write.csv(restable,file=paste0("data/SDM/",paste(sp,collapse="_"),"_summary.csv"),row.names=F)

#colnames(coef)[grepl("X",colnames(coef))]=paste0("Q",sub("X","",colnames(coef)[grepl("X",colnames(coef))]))
#res$param=sub("beta[.]","",rownames(res))
#0resl=melt(res,id.vars=c("model","suitability","param"))
#levels(resl$variable)=c(0.025,0.25,0.5,0.75,0.975)

library(grid) # needed for arrow function
library(ggplot2)

## Make a plot to explore the data


gplot(senv) + geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colours=c('white','blue','red'),na.value="transparent") +
  coord_equal()

## compare predictions
p1=gplot(predr,maxpixels=5e5) + geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colours=c('white','blue','red'),na.value="transparent") +
#  geom_point(data=spd_sp[spd_sp$presence==0,]@data,aes(x=lon,y=lat),pch="-",col="grey",cex=1)+
  geom_point(data=spd_sp[spd_sp$presence==1,]@data,aes(x=lon,y=lat),pch="+",cex=1)+
  coord_equal()+xlab("Longitude")+ylab("Latitude")+
  labs(fill = "p(presence)")+ 
  theme(legend.position="bottom")

p2=ggplot(aucdat, aes(factor(presences>0,labels=c("Absent","Present")), pred))+ 
  geom_boxplot()+facet_wrap(~model)+
  xlab("Observed")+ylab("p(presence)")+
  coord_flip()


## plot the correlogram
p3=ggplot(ac, aes(x=dist2, y=mean, group=model))+
  geom_line(aes(colour=model))+
  geom_pointrange(aes(ymax = max, ymin=min, group=model,colour=model))+
  geom_line(col="black")+
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  ylab("Autocorrelation")+xlab("Distance (km)")+
  theme(legend.position=c(.25, .25))


## make the pdf
pdf(file=paste0("manuscript/figures/SDM_",paste(sp,collapse="_"),".pdf"),width=5,height=10)

grid.newpage()

vp1 <- viewport(width = 1, height = 0.3, x=.5,y=.15)
vp2 <- viewport(width = 1, height = .5, x = .5, y = 0.55)
vp3 <- viewport(width = .95, height = .2, x = .5, y = 0.9)

print(p3, vp = vp1)
print(p1, vp = vp2)
print(p2, vp = vp3)

dev.off()


## Fit a simple GLM to the data
obs=as.matrix(data[,c("presences","trials")])
m_glm=lapply(1:nrow(mods),function(i) glm(paste0("obs",mods$formula[i]),data=data,family="binomial"))

screenreg(m_glm,custom.model.names=as.character(mods$name),format="markdown")

#= Parameter estimates
summary(res[[1]])
xyplot(mod1$mcmc)
#= Predictions
summary(mod1$theta.latent)
summary(mod1$theta.pred)
plot(theta,mod1$theta.pred)

### Check Convergence
xyplot(mod1$mcmc, main="Beta",strip=strip.custom(factor.levels=c("intercept",vars)))
gelman.diag(mod1$mcmc,confidence = 0.95,autoburnin=F,multivariate=T)


