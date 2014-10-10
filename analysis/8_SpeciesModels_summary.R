# Load the libraries and set working directory
source("analysis/setup.R")

library(dismo)
library(grid) # needed for arrow function
library(ggplot2)


#### load summary tables
#ds=do.call(rbind.data.frame,lapply(list.files("output/sdm/",pattern="summary.csv",full=T,recursive=T)[2:3],read.csv))
#ds$species=sub("_"," ",ds$species)
#print(xtable(ds,caption="Evaluation of distribution models using interpolated precipitation or cloud product"), 
#      include.rownames = FALSE,"html",file="manuscript/modeltable.html",digits=2)


## extract environmental data
finputs=list.files("output/sdm/",pattern="modelinput.*.nc$",full=T,recursive=T)
env=hSDM.ncExtract(files=finputs,what="envdata",unscale=T)
idata=filter(hSDM.ncExtract(files=finputs,what="spdata"),trials>0)
idata$value=ifelse(idata$presences>1,1,0)
idata$variable="observations"

## rescale CLD data to 0-1
env[,grepl("CLD",colnames(env))]=env[,grepl("CLD",colnames(env))]/100

## extract model results and evaluation
fresults=list.files("output/sdm/",pattern="modeloutput.*.nc",full=T,recursive=T)

eval=hSDM.ncExtract(files=fresults,what="eval")
coef=hSDM.ncExtract(files=fresults,what="coef")
ac=hSDM.ncExtract(files=fresults,what="autocor")
pred=hSDM.ncExtract(files=fresults,what="predictions")
pred$modelname2=ifelse(pred$modelname=="Cloud","Cloud Model","Precipitation Model")


## combine to a single dataframe
pd=rbind.data.frame(melt(env,id.vars=c("species","x","y","cell"),measure.vars=c("PPTJUL","CLDJUL")),
                    select(pred,c(species,x,y,cell,variable=modelname2,value=pred)),
                    select(idata,c(species,x,y,cell,variable,value)))

## identify regions for each species
pd=mutate(pd,reg=ifelse(species=="Protea_cynaroides"&x>=20&x<=25&y>=-34.8&y<=-33.3|
                          species%in%c("Rupicola_peruvianus","Lepidocolaptes_lacrymiger")&x>=-78&x<=-73&y>=4&y<=5.5,1,0))

#                    left_join(x=pd1,y=select(pred,c(species,cell,modelname,pred)),by=c("species","cell"))

## Graphical Output
pdf(file=paste0("manuscript/figures/SDM_",paste(sp,collapse="_"),".pdf"),width=11,height=7)


library(grid) # needed for arrow function
library(gridExtra)
library(ggplot2)
library(scales)
## Make a plot to explore the data
#fcoast=fortify(crop(coast,ereg))

## characterize species data


penv1=ggplot(filter(pd2,species=="Protea_cynaroides"),aes(x=x,y=y,fill=value)) + geom_tile() +
              facet_grid(variable~species) +
              scale_fill_gradientn(colours=c('white','blue','red'),na.value="transparent")+
              coord_equal(ratio=1.2)+ theme(legend.position="right")+
              ylab("")+xlab("")

penv2=ggplot(filter(pd2,species=="Rupicola_peruvianus"),aes(x=x,y=y,fill=value)) + geom_tile() +
  facet_grid(variable~species) +
  scale_fill_gradientn(colours=c('white','blue','red'),na.value="transparent")+
  coord_equal(ratio=1.2)+ theme(legend.position="right")+
  ylab("")+xlab("")

ppred1=
  ggplot(filter(pred2,species=="Protea_cynaroides"),aes(x=x,y=y,fill=pred)) + geom_tile() +
  facet_grid(~modelname) +
  scale_fill_gradientn(colours=c('white','blue','red'),na.value="transparent")+
  coord_equal(ratio=1.2)+ theme(legend.position="right")+
  ylab("")+xlab("")+
  geom_point(data=filter(idata,trials>1&presences==0&species=="Protea_cynaroides"&x>=20&x<=25&y>=-34.8&y<=-33.3),
             aes(x=x,y=y,fill=1),pch=16,col="green",cex=.1,lwd=2)+
  geom_point(data=filter(idata,trials>1&presences==1&species=="Protea_cynaroides"&x>=20&x<=25&y>=-34.8&y<=-33.3),
             aes(x=x,y=y,fill=1),pch=3,cex=1,lwd=3)+
    labs(fill = "p(presence)")  

ppred2=ggplot(filter(pred2,species=="Rupicola_peruvianus"),aes(x=x,y=y,fill=pred)) + geom_tile() +
  facet_grid(species~modelname) +
  scale_fill_gradientn(colours=c('white','blue','red'),na.value="transparent")+
  coord_equal(ratio=1.2)+ theme(legend.position="right")+
  ylab("")+xlab("")+
  geom_point(data=filter(idata,trials>1&presences==0&species=="Rupicola_peruvianus"&x>=-78&x<=-73&y>=4&y<=5.5),
             aes(x=x,y=y,fill=1),pch=16,col="green",cex=1,lwd=2)+
  geom_point(data=filter(idata,trials>1&presences==1&species=="Rupicola_peruvianus"&x>=-78&x<=-73&y>=4&y<=5.5),
             aes(x=x,y=y,fill=1),pch=3,cex=3,lwd=3)+
  labs(fill = "p(presence)")  


ppred2=ggplot(filter(pred,species=="Rupicola_peruvianus"),aes(x=x,y=y,fill=pred)) + geom_tile() +
  facet_grid(species~modelname) +
  scale_fill_gradientn(colours=c('white','blue','red'),na.value="transparent")+
  coord_equal(ratio=1.2)+ theme(legend.position="right")+
  ylab("")+xlab("")+
  geom_point(data=filter(idata,trials>1&presences==0&species=="Rupicola_peruvianus"),
             aes(x=x,y=y,fill=1),pch=16,col="green",cex=1,lwd=2)+
  geom_point(data=filter(idata,trials>1&presences==1&species=="Rupicola_peruvianus"),
             aes(x=x,y=y,fill=1),pch=3,cex=3,lwd=3)+
  labs(fill = "p(presence)")  


### make one monster?
gs=unique(select(pd,species,variable))
gs$variable=factor(gs$variable,levels=c("PPTJUL","CLDJUL","observations","Precipitation Model","Cloud Model"),ordered=T)
gs=gs[order(gs$variable,gs$species),]
gs$col=as.numeric(as.factor(gs$species))
gs$row=as.numeric(as.factor(gs$variable))
gs$variable=as.character(gs$variable)

plts=lapply(1:nrow(gs),function(i){
  ggplot(filter(pd,reg==1,species==gs$species[i]&variable==gs$variable[i]),aes(x=x,y=y,fill=value)) + geom_tile() +
    facet_grid(species~variable) +
    scale_fill_gradientn(colours=c('white','blue','red'),na.value="transparent")+
    coord_equal(ratio=1.2)+ theme(legend.position="right")+
    ylab("")+xlab("")+
    labs(fill = "")+
    theme(text=element_text(size=48))
})

gs

png(file=paste0("manuscript/figures/SDM_summary_reg.png"),width=3000,height=3000,pointsize=48)
do.call("grid.arrange", c(plts, ncol=max(gs$col)))
dev.off()


plts_full=lapply(1:nrow(gs),function(i){
  ggplot(filter(pd,species==gs$species[i]&variable==gs$variable[i]),aes(x=x,y=y,fill=value)) + geom_tile() +
    facet_grid(species~variable) +
    scale_fill_gradientn(colours=c('white','blue','red'),na.value="transparent")+
    coord_equal(ratio=1.2)+ theme(legend.position="right")+
    ylab("")+xlab("")+
    labs(fill = "")+
    theme(text=element_text(size=48))
})

png(file=paste0("manuscript/figures/SDM_summary_full.png"),width=3000,height=3000,pointsize=48)
do.call("grid.arrange", c(plts_full, ncol=max(gs$col)))
dev.off()


## compare predictions
p1=
  gplot(pred,aes(x=x,y=y,fill=pred,maxpixels=1e4) + 
          geom_tile() +
          facet_grid(species~modeltype) +
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
p3=ggplot(ac, aes(x=dist2, y=mean, group=modelname))+
  geom_line(aes(colour=modelname))+
          facet_grid(.~species)+
  geom_pointrange(aes(ymax = max, ymin=min, group=modelname,colour=modelname))+
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

dev.off()
