# Load the libraries and set working directory
source("analysis/setup.R")

registerDoMC(12)
rasterOptions(datatype="FLT4S")
  
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

data=cbind.data.frame(values(stack(senv,rprot)),coordinates(senv),cell=cellFromXY(senv,xy=coordinates(senv)))
data=as.data.frame(data[!is.na(data[,"presences"]),])
data$trials=as.integer(data$trials)
data$presences=as.integer(data$presences)
data=na.omit(data)

## create 'fitting' dataset where there are observations
fdata=data[data$trials>0,]



### Save model input data
save(env,senv,spr,data,fdata,file=paste0(outputdir,"modelinput.Rdata"))



##############################################################
load(paste0(outputdir,"modelinput.Rdata"))

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


pred$cell=cellFromXY(senv,xy=pred[,c("x","y")])
predr=stack(lapply(unique(pred$model),function(m) rasterFromXYZ(xyz=pred[pred$model==m,c("x","y","pred")])))
names(predr)=unique(pred$model)
projection(predr)='+proj=longlat'
writeRaster(predr*100,file=paste0(outputdir,sp2,".tif"),overwrite=T)

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

##### Autocorrelation
faspac=paste0(outputdir,sp2,"_autocorrelation_FFT.csv")
  ac1=acor_table(predr[[1]],verbose=T)
  ac1$model=names(predr)[1]
  ac2=acor_table(predr[[2]],verbose=T)
  ac2$model=names(predr)[2]
  ac=rbind(ac1,ac2)
write.csv(ac,file=faspac,row.names=F)
}
ac=read.csv(faspac)

#spacs=spac%.%group_by(model)%.%filter(dist2<50)%.%
#  summarize(SCor_mean_50km=median(mean),SCor_sd_50km=diff(range(quantile(mean,c(0.025,0.975)))))

# http://www.bayesian-inference.com/modelfit
# http://voteview.com/DIC-slides.pdf
dic=coef[grepl("Deviance",coef$param),]
dic$pv=(dic$SD^2)/2
dic$dic=dic$Mean+dic$pv
dic

## AUC
aucdat=merge(fdata[,c("presences","trials","cell")],pred,by=c("cell"))
#aucdat$model=factor(aucdat$model,levels=mods$model,ordered=T)

feval=function(x){
  e=evaluate(p=x$pred[x$presences>0],a=x$pred[x$presences==0])
  data.frame(nPresence=e@np,nTrials=e@na,auc=e@auc,cor=e@cor)
}

restable=do.call(rbind.data.frame,by(aucdat,INDICES=aucdat$modelname,FUN=feval))
restable$species=sp2
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


#colnames(coef)[grepl("X",colnames(coef))]=paste0("Q",sub("X","",colnames(coef)[grepl("X",colnames(coef))]))
#res$param=sub("beta[.]","",rownames(res))
#0resl=melt(res,id.vars=c("model","suitability","param"))
#levels(resl$variable)=c(0.025,0.25,0.5,0.75,0.975)

library(grid) # needed for arrow function
library(ggplot2)

## Make a plot to explore the data
pdf(file="manuscript/figures/SDM.pdf",width=11,height=7)

gplot(senv) + geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colours=c('white','blue','red'),na.value="transparent") +
  coord_equal()

ggplot(coef[!coef$param%in%c("Vrho","Deviance"),], aes(ymin = X2.5.,ymax=X97.5.,y=X50., x = param,colour=modelname)) + 
  geom_hline(yintercept=0)+
  geom_point(position = position_dodge(.5))+
  geom_linerange(position = position_dodge(.5),lwd=1)+
  ylab("Value")+xlab("Parameter")+
  coord_flip()

ggplot(ac, aes(x=dist2, y=mean, group=model))+
  geom_line(aes(colour=model))+
  geom_pointrange(aes(ymax = max, ymin=min, group=model,colour=model))+
  geom_line(col="black")+
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  ylab("Autocorrelation")+xlab("Distance (km)")+
  theme(legend.position=c(.25, .25))

ggplot(pred) + geom_tile(aes(x=x,y=y,fill = pred)) +
  facet_wrap(~modelname,ncol=1) +
  scale_fill_gradientn(colours=c('white','blue','red')) +
  geom_point(data=prot[prot$pro==sp,]@data,aes(x=londd,y=latdd),pch="+")+
  xlim(c(17.8,25.2))+ylim(-35,-32)+
  coord_equal()

#ggplot(rho) + geom_tile(aes(x=x,y=y,fill = rho)) +
#  facet_wrap(~suitability,ncol=1) +
#  scale_fill_gradientn(colours=c('white','blue','red')) +
#  coord_equal()

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
  
p2=gplot(env[["ALT"]],maxpixels=5e6) + geom_tile(aes(fill = value)) +
  facet_grid(~variable,labeller=function(...) return("Protea cynaroides occurrences")) +
    scale_fill_gradientn(colours=c('darkgreen','yellow','red'),na.value="transparent",
                       name="Elevation (m)") +
  coord_equal()+
  geom_point(data=prot[prot$pro!=sp,]@data,aes(x=londd,y=latdd),col="grey",alpha=.2,pch=1,cex=.5)+
  geom_point(data=prot[prot$pro==sp,]@data,aes(x=londd,y=latdd),pch="+",cex=2)+
  scale_colour_manual(name = 'Occurrence', 
                      values =c('black'='black','grey'='grey'),
                      labels = c('Presence','Absence'))+
  xlim(c(20,25))+ylim(-34.2,-33.5)+
  xlab(label="")+
  theme(panel.background = element_rect(fill='transparent'),legend.key.width = unit(1.5, "cm"))+
  ylab("")


grid.newpage()
vp1 <- viewport(width = 1, height = 0.67, x=.5,y=.33)
vp2 <- viewport(width = 1, height = .33, x = .5, y = 0.72)
print(p1, vp = vp1)
print(p2, vp = vp2)
dev.off()

gplot(rprot[["presences"]]) + geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colours=c('white','blue','red'),na.value="transparent") +
  coord_equal()
