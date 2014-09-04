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
#taxon=taxon%.%filter(GENUS_NAME==sp[1])
taxon=taxon%.%filter(FAMILY_NAME==fam)

idcols=c("SAMPLING_EVENT_ID","LATITUDE","LONGITUDE","YEAR","MONTH","DAY","TIME","COUNTRY","STATE_PROVINCE","COUNT_TYPE","EFFORT_HRS","EFFORT_DISTANCE_KM",
         "EFFORT_AREA_HA","OBSERVER_ID","NUMBER_OBSERVERS","GROUP_ID","PRIMARY_CHECKLIST_FLAG")
spcols=as.character(taxon$SCI_NAME)

library(data.table)
library(ff)
library(sqldf)

f=list.files(ebirddir,pattern="checklists.csv",recursive=T,full=T)[2]

d=foreach(f=list.files(ebirddir,pattern="checklists.csv",recursive=T,full=T),.combine=rbind) %dopar% {
  ## read in data
#  f=list.files(ebirddir,pattern="checklists.csv",recursive=T,full=T)[10]
#  d=fread( f,header = T, sep = ',',select=c(idcols,spcols),na.strings="?",showProgress=T,verbose=T)
#  d=read.csv.ffdf(file=f)
  
  ## get colnames in table
  hdr=colnames(read.csv(f,nrows=1))
## which colnames are in our list
  tsp=grep(paste(c(idcols,spcols),collapse="|"),hdr,value=T)
   d=system.time(read.csv.sql(f, sql = paste0("select \"",paste(tsp,collapse="\",\""),"\" FROM file WHERE ",
                                 "LATITUDE  BETWEEN ",bbox(reg)["y","min"]," AND ",bbox(reg)["y","max"]," AND ",
                                 "LONGITUDE BETWEEN ",bbox(reg)["x","min"]," AND ",bbox(reg)["x","max"]," AND ",
                                 "PRIMARY_CHECKLIST_FLAG=1",";")))
  colnames(d)=gsub("\"","",colnames(d))
  ## loop through species and clean up "X" flags
  for(i in tsp){
    x=d[,i]
    x[x=="X"]="1"
    x=as.numeric(x)
    x[x>0]=1
    d[,i]=x
    #  print(i)
  }
  ## calculate number of trials and presences
  d$trials=rowSums(d[,tsp])>0
  d=d%.%filter(trials>0)
  d$presence=d[,paste(sp,collapse="_")]
  ## subset to columns of interest
  d=d[,c("LATITUDE","LONGITUDE","YEAR","MONTH","presence","trials")]
  print(f)
  return(d)
}

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

sp="PRCYNA"
rprot=fprot(prot,sp=sp)

senv=scale(env)

#splom(senv)

data=cbind.data.frame(values(stack(senv,rprot)),coordinates(senv),cell=cellFromXY(senv,xy=coordinates(senv)))
data=as.data.frame(data[!is.na(data[,"presences"]),])
data$trials=as.integer(data$trials)
data$presences=as.integer(data$presences)
data=na.omit(data)

## adjacency matrix
neighbors.mat <- adjacent(senv, cells=1:ncell(senv), directions=8, pairs=TRUE, sorted=TRUE)
#neighbors.mat=neighbors.mat[neighbors.mat[,1]%in%data$cell,]
n.neighbors <- as.data.frame(table(as.factor(neighbors.mat[,1])))[,2]
adj <- neighbors.mat[,2]

## create 'fitting' dataset where there are observations
fdata=data[data$trials>0,]

nchains=3

mods=data.frame(
  model=c("m1","m2","m3"),
  formula=c("~tmean1+tmean7+prec1+prec7+alt",
                "~tmean1+tmean7+c_mm_01+c_mm_07+alt",
                "~tmean1+tmean7+prec1+prec7+c_mm_01+c_mm_07+alt"),
  name=c("Temperature & Precipitation",
         "Temperature & Cloud",
         "Temperature, Precipitation, & Cloud"))

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
#      spatial.entity=fdata$cell,
#      n.neighbors=n.neighbors,
#      neighbors=adj,
#      spatial.entity.pred=data$cell,
      burnin=1000, mcmc=5000, thin=5,
      beta.start=0,
      suitability.pred=data,
#      Vrho.start=20,
#      priorVrho="1/Gamma",
#      shape=1, rate=1,
#      save.rho=0, 
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

coef=do.call(rbind,lapply(res,function(x) x$coef))
pred=do.call(rbind,lapply(res,function(x) x$pred))
rho=do.call(rbind,lapply(res,function(x) x$rho))

pred=pred[pred$model%in%mods$model[1:2],]

# http://www.bayesian-inference.com/modelfit
# http://voteview.com/DIC-slides.pdf
dic=coef[grepl("Deviance",coef$param),]
dic$pv=(dic$SD^2)/2
dic$dic=dic$Mean+dic$pv
dic

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

ggplot(coef[coef$param=="Deviance",], aes(ymin = X2.5.,ymax=X97.5.,y=X50., x = param,colour=modelname)) + 
  geom_point(position = position_dodge(.5))+
  geom_linerange(position = position_dodge(.5),lwd=1)+
  ylab("Value")+xlab("Parameter")+
  coord_flip()

ggplot(coef[!coef$param%in%c("Vrho","Deviance"),], aes(ymin = X2.5.,ymax=X97.5.,y=X50., x = param,colour=modelname)) + 
  geom_hline(yintercept=0)+
  geom_point(position = position_dodge(.5))+
  geom_linerange(position = position_dodge(.5),lwd=1)+
  ylab("Value")+xlab("Parameter")+
  coord_flip()

#ggplot(coef[coef$param%in%c("Vrho"),], aes(ymin = X2.5.,ymax=X97.5.,y=X50., x = param,colour=modelname)) + 
#  geom_hline(yintercept=0)+
#  geom_point(position = position_dodge(.5))+
#  geom_linerange(position = position_dodge(.5),lwd=1)+
#  ylab("Value")+xlab("Parameter")+
#  coord_flip()

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
  
p2=gplot(env[["alt"]],maxpixels=5e6) + geom_tile(aes(fill = value)) +
  facet_grid(~variable,labeller=function(...) return("Protea cynaroides occurrence")) +
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



densityplot(mod1$mcmc)#, main="Posterior Distributions",
            strip=strip.custom(factor.levels=c("intercept",vars)),
            scales=list(relation="same"),layout=c(1,7))+
  layer(panel.abline(v=0))
HPDinterval(ps.beta[[1]], prob = 0.95)

## Predict model to the grid

```{r,predictmodel}
## First subset area to speed up predictions
pext=extent(c(-50,-48,-26.5,-24))
penv=crop(senv,pext)

## if you want to make predictions for the full grid, run this line:
#penv=senv

## Calculate posterior estimates of p(occurrence) for each cell
## This extracts the posterior coefficients, performs the regression, 
## calculates the quantiles, and takes the inverse logit to get p(occurrence)

## niter will use a reduced number of posterior samples to generate the summaries
pred=calc(penv,function(x,niter=30) {
  mu1=apply(apply(ps.beta[[1]][1:niter,],1,function(y) y*c(1,x)),2,sum,na.rm=T)
  mu2=quantile(mu1,c(0.025,0.5,0.975),na.rm=T)  
  p=1/(1+exp(-mu2))
  return(p)
})
names(pred)=c("Lower_CI_2.5","Median","Upper_CI_97.5")
## Write out the predictions
writeRaster(pred,file="Prediction.tif",overwrite=T)
```
