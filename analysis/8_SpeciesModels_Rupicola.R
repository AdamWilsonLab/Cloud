# Load the libraries and set working directory
source("analysis/setup.R")
install_github(repo="Rdatatable/data.table")
library(data.table)
library(AUC)
library(dismo)
library(redshift)
library(rasterVis)

#library(devtools)
#install_github("adammwilson/rasterAutocorr")
library(rasterAutocorr)

ncores=12
registerDoMC(ncores)
rasterOptions(datatype="FLT4S")

#sp=c("Lepidocolaptes","lacrymiger")
#fam="Furnariidae (Ovenbirds and Woodcreepers)"

## cock of the rock
sp=c("Rupicola","peruvianus")
fam="Cotingidae (Cotingas)"

sp2=paste0(sp,collapse="_")

### Load eBird data for 
## ebird data downloaded and unzipped from http://ebirddata.ornith.cornell.edu/downloads/erd/ebird_all_species/erd_western_hemisphere_data_grouped_by_year_v5.0.tar.gz
ebirddir="/mnt/data/jetzlab/Data/specdist/global/ebird/"

## select species for presences and 'absences'
taxon=read.csv(paste0(ebirddir,"/doc/taxonomy.csv"),na.strings="?")
taxon$TAXON_ORDER2=round(taxon$TAXON_ORDER)
## some duplicate species names marked with "/"
## Many SPECIES_NAME and GENUS_NAME are null but SCI_NAME is not
taxon$gensp=apply(do.call(rbind,lapply(strsplit(as.character(taxon$SCI_NAME),split="_|/"),function(x) x[1:2])),1,paste,collapse="_")
taxon=taxon%.%filter(FAMILY_NAME==fam)

## Select list (or single) taxon ids for use as 'presence'
sptaxon=taxon$TAXON_ORDER2[taxon$gensp%in%sp2]
## Select list of taxon ids for use as 'non-detection/absence'
nulltaxon=taxon$TAXON_ORDER2[taxon$TAXON_ORDER2!=sptaxon]

## create output folder
outputdir=paste0("output/sdm/",sp2,"/")
if(!file.exists(outputdir)) dir.create(outputdir,recursive=T)

## load region boundary
download.file(paste0("http://mol.cartodb.com/api/v2/sql?",
                  "q=SELECT%20ST_TRANSFORM(the_geom_webmercator,4326)%20as%20the_geom,%20seasonality%20FROM%20",
                  "get_tile('jetz','range','",
                  sp[1],"%20",sp[2],
                  "','jetz_maps')&format=shp&filename=expert"),
              destfile=paste0(outputdir,"expert.zip"))
unzip(paste0(outputdir,"expert.zip"),exdir=outputdir)



reg=readOGR(outputdir,"expert")
ereg=extent(reg)

## adjust bbox if desired
ereg@xmin=-81.4

## get species data
t1=system.time(spd<<-getebird(con=conn,sptaxon=sptaxon,nulltaxon=nulltaxon,region=reg))
coordinates(spd)=c("longitude","latitude")
projection(spd)="+proj=longlat +datum=WGS84 +ellps=WGS84"
spd@data[,c("lon","lat")]=coordinates(spd)  

## print a little summary

writeLines(paste(sp2," has ",nrow(spd),"rows"))

#########  Environmental Data
cf=stack(c("data/MCD09/MCD09_mean_01.tif","data/MCD09/MCD09_mean_07.tif","data/MCD09_deriv/intra.tif"))
names(cf)=c("CLDJAN","CLDJUL","CLDSEAS")
NAvalue(cf)=0
cf=crop(cf,ereg)
gain(cf)=.01

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

## rasterize species data to environmental grid
presences=rasterize(spd[spd$presence==1,],cf,fun="count",field="presence",background=0)
trials=rasterize(spd[spd$presence==0,],cf,fun="count",field="presence",background=0)
spr=stack(presences,trials)
names(spr)=c("presences","trials")

senv=scale(env)


data=cbind.data.frame(values(stack(senv,spr)),coordinates(senv),cell=cellFromXY(senv,xy=coordinates(senv)))
data=as.data.frame(data[!is.na(data[,"presences"]),])
data$trials=as.integer(data$trials)
data$presences=as.integer(data$presences)
data=na.omit(data)

## create 'fitting' dataset where there are observations
fdata=data[data$trials>0,]
fdata=data[data$trials>0&data$trials>=data$presences,]

## due to opportunistic observations, there are a few sites with more presences than trials, update those records so presences=trials
fdata$trials[fdata$presences>fdata$trials]=fdata$presences[fdata$presences>fdata$trials]

### Save model input data
save(env,senv,spr,data,fdata,file=paste0(outputdir,"modelinput.Rdata"))

#####################################################################################
#####################################################################################
load(paste0(outputdir,"modelinput.Rdata")))

nchains=3

mods=data.frame(
  model=c("m1","m2"),
  formula=c("~PPTJAN+PPTJUL+PPTSEAS+MAT",
            "~CLDJAN+CLDJUL+CLDSEAS+MAT"),
  name=c("Precipitation",
         "Cloud"))

## file to save output
fres=paste0(outputdir,"modeloutput.Rdata")

registerDoMC(ncores)

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
return(list(coef=coef,pred=pred,rho=rho))  
}
save(res,file=fres)
}

#################
load(fres)

coef=do.call(rbind,lapply(res,function(x) x$coef))
pred=do.call(rbind,lapply(res,function(x) x$pred))

#pred=pred[pred$model%in%mods$model[c(1,3)],]
pred$cell=cellFromXY(senv,xy=pred[,c("x","y")])

predr=stack(lapply(unique(pred$model),function(m) rasterFromXYZ(xyz=pred[pred$model==m,c("x","y","pred")])))
names(predr)=unique(pred$model)
projection(predr)='+proj=longlat'
writeRaster(predr*100,file=paste0(outputdir,sp2,".tif"),overwrite=T)


##### Autocorrelation
ac1=acor_table(predr[[1]],verbose=T)
ac1$model=names(predr)[1]
ac2=acor_table(predr[[2]],verbose=T)
ac2$model=names(predr)[2]

ac=rbind.data.frame(ac1,ac2)

## observed data
mask=ifelse(trials>0,1,0)
NAvalue(mask)=0
pres1=mask(presences,mask)
ac_obs=acor_table(pres1,verbose=T)
ac_obs$model="Observed"


## Global Spatial Autocorrelation
fspac=paste0(outputdir,sp2,".csv")
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

# http://www.bayesian-inference.com/modelfit
# http://voteview.com/DIC-slides.pdf
dic=coef[grepl("Deviance",coef$param),]
dic$pv=(dic$SD^2)/2
dic$dic=dic$Mean+dic$pv

## AUC
aucdat=merge(fdata[,c("presences","trials","cell")],pred,by=c("cell"))
#aucdat$model=factor(aucdat$model,levels=mods$model,ordered=T)

feval=function(x){
  e=evaluate(p=x$pred[x$presences>0],a=x$pred[x$presences==0])
  data.frame(nPresence=e@np,nTrials=e@na,auc=e@auc,cor=e@cor)
}

restable=do.call(rbind.data.frame,by(aucdat,INDICES=aucdat$modelname,FUN=feval))
restable$species=paste(sp,collapse=" ")
restable$model=dic$model[match(rownames(restable),dic$modelname)]
restable$Deviance=dic$Mean[match(restable$model,dic$model)]
restable$Pv=dic$pv[match(restable$model,dic$model)]
restable$DIC=dic$dic[match(restable$model,dic$model)]
##
restable[,c("MoransI","GearyC")]=spac[match(restable$model,spac$model),c("MoransI","GearyC")]

restable=restable[order(restable$model),
                  c("species","model","nPresence","nTrials","Deviance","Pv","DIC","auc","cor","MoransI","GearyC")]

write.csv(restable,file=paste0(outputdir,sp2,"_summary.csv"),row.names=F)

print(xtable(restable,caption="Evaluation of distribution models using interpolated precipitation or cloud product"), 
      include.rownames = FALSE,"html",file=paste0(outputdir,sp2,"_summary.html"))


 ## add elevation gradient

## Correlogram of points
#library(ncf)
#s=sample(1:nrow(fdata),100)
#pcg=correlog(x=fdata$x[s],y=fdata$y[s],z=fdata$presence[s]>0,latlon=T,increment=50)

#library(spdep)
#ds=c(seq(0,100,by=5))
#joincount.mc(HICRIME, nb2listw(COL.nb, style="B"), nsim=99)
#d1=dnearneigh(as.matrix(fdata[,c("x","y")]),d1=0,d2=100,longlat=T)
#plot(pcg);abline(h=0,col="red")

#colnames(coef)[grepl("X",colnames(coef))]=paste0("Q",sub("X","",colnames(coef)[grepl("X",colnames(coef))]))
#res$param=sub("beta[.]","",rownames(res))
#0resl=melt(res,id.vars=c("model","suitability","param"))
#levels(resl$variable)=c(0.025,0.25,0.5,0.75,0.975)
library(spdep)
set.ZeroPolicyOption(T)  

nb=dnearneigh(spd,0,5,longlat=T)
plot(nb)
nblist=nb2listw(nb,style="B")
#tdrop=do.call(c,lapply(nblist$weights,function(x) length(x)==0))
#nblist=nb2listw(nb,style="B",zero.policy=TRUE)
joincount.test(as.factor(spd$presence>0),nblist)

library(grid) # needed for arrow function
library(ggplot2)
library(scales)
## Make a plot to explore the data
#fcoast=fortify(crop(coast,ereg))

penv=gplot(senv) + geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colours=c('white','blue','red'),na.value="transparent") +
  coord_equal()+ theme(legend.position="bottom")

## compare predictions
p1=
  gplot(predr,maxpixels=1e4) + geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
   scale_fill_gradientn(colours=c('white','blue','red'),values = c(0,.6,.75,1), 
                        rescaler = function(x, ...) x,
                       na.value="transparent") +
  coord_equal()+xlab("Longitude")+ylab("Latitude")+
  geom_point(data=spd[spd$presence==0,]@data,aes(x=lon,y=lat),pch=1,col=grey(.8),cex=1)+
  geom_point(data=spd[spd$presence==1,]@data,aes(x=lon,y=lat),pch=3,cex=2)+
#  geom_line(data=fcoast,mapping=aes(long,lat,group=group)) +
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

p4=  gplot(predr,maxpixels=1e7) + geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colours=c('white','blue','red'),values = c(0,.6,.75,1), 
                       rescaler = function(x, ...) x,
                       na.value="white") +
  coord_equal()+xlab("Longitude")+ylab("Latitude")+
  geom_point(data=spd[spd$presence==0,]@data,aes(x=lon,y=lat),pch=1,col=grey(.8),cex=1)+
  geom_point(data=spd[spd$presence==1,]@data,aes(x=lon,y=lat),pch=3,cex=2)+
  labs(fill = "p(presence)")+ 
  theme(legend.position="bottom")+
  ylim(c(-1,7.5))+xlim(c(-79.5,-72))

## make the pdf
pdf(file=paste0(outputdir,"env_",sp2,".pdf"),width=5,height=10)
penv
dev.off()


pdf(file=paste0(outputdir,"SDM_",sp2,".pdf"),width=5,height=10)
grid.newpage()

vp1 <- viewport(width = 1, height = 0.3, x=.5,y=.15)
vp2 <- viewport(width = 1, height = .5, x = .5, y = 0.55)
vp3 <- viewport(width = .95, height = .2, x = .5, y = 0.9)

print(p3, vp = vp1)
print(p4, vp = vp2)
print(p2, vp = vp3)

dev.off()



ggplot(coef[coef$param!="Deviance",], aes(x = Mean, y = param,colour=model))+  
  geom_segment(aes(xend = X97.5.,yend=param),lwd=2)

### Check Convergence
#xyplot(mod1$mcmc, main="Beta",strip=strip.custom(factor.levels=c("intercept",vars)))
#gelman.diag(mod1$mcmc,confidence = 0.95,autoburnin=F,multivariate=T)


