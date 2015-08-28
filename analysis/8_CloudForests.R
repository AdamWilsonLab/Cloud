
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
env=crop(stack(bprods),extent(c(-160,160,-23.43726,23.43726))) ## all tropics

wc=crop(stack(c(
  "/mnt/data/jetzlab/Data/environ/global/worldclim/alt.bil",
  "/mnt/data/jetzlab/Data/environ/global/worldclim/bio_1.bil",
  "/mnt/data/jetzlab/Data/environ/global/worldclim/bio_12.bil",
  "/mnt/data/jetzlab/Data/environ/global/worldclim/bio_15.bil"
)),env)
names(wc)=c("elev","MAT","MAP","PSeas")

env2=stack(env,wc)


## load land map
land=map(interior=F,fill=T,xlim=bbox(env2)[1,],ylim=bbox(env2)[2,],plot=F)
land=map2SpatialPolygons(land,IDs=1:length(land$names))
## mask ocean
env3=mask(env2,land,file="data/out/CloudEnv.tif",options=c("COMPRESS=LZW"),overwrite=T)
env4=scale(env3)
region=raster(env[[1]])
values(region)=as.factor(ifelse(coordinates(region)[,"x"]<(-29.5),"Americas",ifelse(coordinates(region)[,"x"]<63.4,"Africa","Asia Pacific")))
env4=stack(env4,region)
 writeRaster(env4,file="data/out/CloudEnv_scaled.tif",options=c("COMPRESS=LZW"),overwrite=T)
}

env=stack("data/out/CloudEnv_scaled.tif")
names(env)=sub("[.]tif","",c(basename(bprods),names(wc),"region"))
levels(env[["region"]])=data.frame(ID=1:3,code=c("Africa","Americas","Asia Pacific"))

pres=na.omit(cbind.data.frame(cf=1,coordinates(cfp),raster::extract(env,cfp,df=T,ID=F)))
colnames(pres)[2:3]=c("x","y")

ns=10000
abs=cbind.data.frame(cf=0,raster::sampleRandom(env,ns,xy=T))
abs$ID=1:nrow(abs)

data=bind_rows(pres,abs)

## add flag for region
data$regionname=as.factor(ifelse(data$x<(-29.5),"Americas",ifelse(data$x<63.4,"Africa","Asia Pacific")))

save(data,file="data/out/CloudForestPoints.Rdata")


################################################
load("data/out/CloudForestPoints.Rdata")
env=stack("data/out/CloudEnv_scaled.tif")
names(env)=sub("[.]tif","",c(basename(bprods),names(wc),"region"))
levels(env[["region"]])=data.frame(ID=1:3,code=c("Africa","Americas","Asia Pacific"))

#data$region=factor(data$region,labels=c("Africa","Americas","AsiaPacific"))
data$region=as.factor(data$region)

models=rbind.data.frame(
  c("Interpolated","cf ~ elev+I(elev^2)+MAT+I(MAT^2)+MAP+I(MAP^2)+PSeas+region*MAP"),
  c("Cloud","cf ~ elev+I(elev^2)+inter+intra+meanannual*region"))
  
colnames(models)=c("name","formula")

mods=foreach(f=models$formula[!grepl("step",models$name)]) %dopar%{
  glm(as.formula(as.character(f)), family=binomial,data=data,weights=1E3^(1-cf))
}

##mods[["stepworldclim"]]=step(mods[[grep("^worldclim$",models$name,fixed=F)]],scope="cf~.^2")
#mods[["stepcloud"]]=step(mods[[grep("^cloud$",models$name,fixed=F)]],scope="cf~.^2")
#mods[["stepall"]]=step(mods[[grep("^all$",models$name)]],scope="cf~.^2")

names(mods)=models$name

screenreg(mods,bold = 0.05,digits=3,single.row=T,
          reorder.coef=c(1:8,13:15,9:12,16:17),
          groups = list("WorldClim" = 4:8, "Cloud Product"=9:11, "Regions" = 12:13,"Interactions"=14:17))

htmlreg(mods,file="output/CloudForest.html",bold = 0.05,digits=2,stars=c(0.001, 0.01, 0.05),caption="Regression Summary",single.row=T,
        reorder.coef=c(1:8,13:15,9:12,16:17))

### compare models
BIC(mods[[1]])-BIC(mods[[2]])

#visreg(mods$step)
#visreg2d(mods$step,"meanannual","intra")
#visreg(mods$step,xvar="meanannual","intra")

beginCluster(10)
mcoptions <- list(preschedule=FALSE, set.seed=FALSE)

ptype="response"

psi=1:nrow(models)

ps=foreach(i=1:nrow(models),.options.multicore=mcoptions,.combine=stack)%dopar%{
  fo=paste0("data/out/CloudForestPrediction_",models$name[i],".tif")
  p1=predict(env,mods[[i]],type=ptype,file=fo,overwrite=T,factors=list(region=c("Africa","Americas","AsiaPacific")))
  raster(fo)
  }


#### Read in predictions

ps=stack(list.files("data/out/",pattern="CloudForestPrediction_Cloud Product.tif",full=T))

#psd=ps[[1]]-ps[[2]]

hcols=function(x,bias=1) {
  #colorRampPalette(c('grey20','grey30','grey40','grey50','steelblue4','steelblue2','gold','red1','red4','red4'),bias=bias)(x) 
  colorRampPalette(c('white','grey90','grey80','grey70','steelblue2','gold','red1','red4','red4'),bias=bias)(x) 
  
}

p_p1=
  gplot(ps,max=1e6)+
  geom_raster(aes(fill=value))+#facet_wrap("variable",ncol=1)+
  scale_fill_gradientn(colours=hcols(1000,bias=.35),trans = "log10",#lim=c(1e-3,1),
                       name="Relative\nOccurrence\nRate\np(x|Y=1)",na.value="transparent")+
  coord_equal(xlim=c(-100,160))+
#  geom_polygon(aes(x=long,y=lat,group=group),
#               data=fortify(land),
#               fill="transparent",col="black",size=.2)+
  geom_point(aes(x = x, y = y), 
             data = data[data$cf==1,],
             col="black",size=.5,shape=1)+
  ylab("")+xlab("")

png("figure/CloudForest.png",width=3000,height=1000,res=300)
print(p_p1)
grid.edit("geom_point.points", grep = TRUE, gp = gpar(lwd = .3)) #make points thinner
dev.off()


## Calculate thresholded area
thresh=function(model) {
  evaluate(model=model,p=model$fitted.values[model$model[,response]==1],a=model$fitted.values[model$model[,response]==0])
  threshold(e1)
}

## calculate area of each cell
area=area(env[[1]])

## get threshold
t1=thresh(mods[["Cloud Product"]])$kappa
## map the threshold
pst=ps>=t1
## get area for each TMCF pixel
psta=pst*area

## calculate total area of tropical cloud forests in km^2
cellStats(pst,sum)

p_p2=
  gplot(pst,max=1e6)+
  geom_raster(aes(fill=value))+#facet_wrap("variable",ncol=1)+
  coord_equal(xlim=c(-100,160))+
#  geom_polygon(aes(x=long,y=lat,group=group),
#               data=fortify(land),
#               fill="transparent",col="black",size=.2)+
  geom_point(aes(x = x, y = y), 
             data = data[data$cf==1,],
             col="black",size=.5,shape=1)+
  ylab("")+xlab("")

png("figure/CloudForest_Thresholded.png",width=3000,height=1000,res=300)
print(p_p2)
grid.edit("geom_point.points", grep = TRUE, gp = gpar(lwd = .2)) #make points thinner
dev.off()


#gplot(psd,max=1e5)+geom_raster(aes(fill=value))+
#  scale_fill_gradient2(low = "red", mid = "white", high = "blue",
#                        midpoint = 0, space = "rgb", na.value = "grey50", guide = "colourbar")+
#  coord_equal(xlim=c(90,160))+
  #geom_polygon(aes(x=long,y=lat,group=group),
  #             data=fortify(range),
  #             fill="transparent",col="black",size=.2)+
#  geom_point(aes(x = x, y = y), 
#             data = data[data$cf==1,],
#             col="black",size=2)

d1=fortify(mods[[2]])

ggplot(d1,aes(y=.hat,x=elev))+geom_line()
ggplot(d1,aes(y=.hat,x=meanannual))+geom_point()+geom_smooth()+scale_y_log10()

