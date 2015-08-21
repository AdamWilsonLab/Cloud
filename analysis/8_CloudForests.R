
##  Cloud forests explanation
# Options: use expert based thresholding and points as validation
# use points to fit logistic for region and validation with additional polygon
# able to monitor, not a one-off map

library(dplyr)
library(AUC)
library(texreg)
library(visreg)
library(raster)

library(doParallel)
registerDoParallel(5)

cfp=readOGR("data/src//cloud_forests/","cloud_forest_points_1997")

bprods=c("data/MCD09_deriv/inter.tif",
         "data/MCD09_deriv/intra.tif",
         "data/MCD09_deriv/seasconc.tif",
         "data/MCD09_deriv/seastheta.tif",
         "data/MCD09_deriv/meanannual.tif",
         "data/MCD09_deriv/mean_1deg_sd.tif")

if(!file.exists("data/out/CloudEnv.tif")){

env=crop(stack(bprods),extent(c(90,171,-14,20)))
env=crop(stack(bprods),extent(c(-160,160,-23.43726,23.43726))) ## global

wc=crop(stack(c(
  "/mnt/data/jetzlab/Data/environ/global/worldclim/alt.bil",
  "/mnt/data/jetzlab/Data/environ/global/worldclim/bio_1.bil",
  "/mnt/data/jetzlab/Data/environ/global/worldclim/bio_12.bil",
  "/mnt/data/jetzlab/Data/environ/global/worldclim/bio_15.bil"
)),env)

env2=stack(env,wc)


## load land map
land=map(interior=F,fill=T,xlim=bbox(env2)[1,],ylim=bbox(env2)[2,],plot=F)
land=map2SpatialPolygons(land,IDs=1:length(land$names))
## mask ocean
env3=mask(env2,land,file="data/out/CloudEnv.tif",options=c("COMPRESS=LZW"),overwrite=T)
env4=scale(env3)
 writeRaster(env4,file="data/out/CloudEnv_scaled.tif",options=c("COMPRESS=LZW"),overwrite=T)
}

env=stack("data/out/CloudEnv_scaled.tif")
names(env)=sub("[.]tif","",c(basename(bprods),"elev"))

pres=na.omit(cbind.data.frame(cf=1,coordinates(cfp),raster::extract(env,cfp,df=T,ID=F)))
colnames(pres)[2:3]=c("x","y")

ns=10000
abs=cbind.data.frame(cf=0,raster::sampleRandom(env,ns,xy=T))
abs$ID=1:nrow(abs)

data=bind_rows(pres,abs)

## add flag for region
data$region=as.factor(ifelse(data$x<(-29.5),"Americas",ifelse(data$x<63.4,"Africa","Asia Pacific")))

save(data,file="data/out/CloudForestPoints.Rdata")


################################################
load("data/out/CloudForestPoints.Rdata")
env=stack("data/out/CloudEnv_scaled.tif")
names(env)=sub("[.]tif","",c(basename(bprods),"elev"))


models=rbind.data.frame(
  c("elevation","cf ~ elev+I(elev^2)"),
  c("cloud+elevation","cf ~ inter+I(inter^2)+elev+I(elev^2)+intra+I(intra^2)+meanannual+I(meanannual^2)"),
  c("all","cf ~ elev+I(elev^2)+inter+I(inter^2)+intra+I(intra^2)+meanannual+I(meanannual^2)+inter*intra*meanannual+region"),
  c("justcloud","cf ~ inter+I(inter^2)+intra+I(intra^2)+meanannual+I(meanannual^2)"),
  c("step","StepwiseSelection"))

colnames(models)=c("name","formula")

mods=foreach(f=models$formula[-nrow(models)]) %dopar%{
  glm(as.formula(as.character(f)), family=binomial,data=data,weights=1E3^(1-cf))
}

mods[[nrow(models)]]=step(mods[[3]])

models$AIC=unlist(lapply(mods,function(x) AIC(x)))
models$AUC=unlist(lapply(mods,function(x) {
  r1=roc(x$fitted.values,as.factor(x$data$cf))
  auc(r1)
}))

names(mods)=models$name

screenreg(mods)

visreg(mods$step)
visreg2d(mods$step,"meanannual","intra")
visreg(mods$step,xvar="meanannual","intra")

#beginCluster(10)
mcoptions <- list(preschedule=FALSE, set.seed=FALSE)

ptype="response"

psi=1:nrow(models)

ps=foreach(i=1:nrow(models),.options.multicore=mcoptions,.combine=stack)%dopar%{

  fo=paste0("data/out/CloudForestPrediction_",models$name[i],".tif")
  p1=predict(env,mods[[i]],type=ptype,file=fo,overwrite=T)
  raster(fo)
  }

ps=stack(list.files("data/out/",pattern="CloudForestPrediction_step.tif",full=T))

psd=ps[[1]]-ps[[2]]

hcols=function(x,bias=1) {
  colorRampPalette(c('grey90','steelblue4','steelblue2','steelblue1','gold','red1','red4'),bias=bias)(x) 
}

p_p1=
  
  gplot(ps,max=5e6)+
  geom_raster(aes(fill=value))+facet_wrap("variable",ncol=1)+
  scale_fill_gradientn(colours=hcols(1000,bias=.2),trans = "log10",#lim=c(1e-3,1),
                       name="Relative Occurrence Rate\np(x|Y=1)",na.value="transparent")+
  coord_equal(xlim=c(-100,160))+
  #coord_equal(xlim=c(90,160))+
#  geom_polygon(aes(x=long,y=lat,group=group),
#               data=fortify(land),
#               fill="transparent",col="black",size=.2)+
  geom_point(aes(x = x, y = y), 
             data = data[data$cf==1,],
             col="black",size=1,solid = FALSE)

gplot(psd,max=1e5)+geom_raster(aes(fill=value))+
  scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                        midpoint = 0, space = "rgb", na.value = "grey50", guide = "colourbar")+
  coord_equal(xlim=c(90,160))+
  #geom_polygon(aes(x=long,y=lat,group=group),
  #             data=fortify(range),
  #             fill="transparent",col="black",size=.2)+
  geom_point(aes(x = x, y = y), 
             data = data[data$cf==1,],
             col="black",size=2)

d1=fortify(mods[[2]])

ggplot(d1,aes(y=.hat,x=elev))+geom_line()
ggplot(d1,aes(y=.hat,x=meanannual))+geom_line()

