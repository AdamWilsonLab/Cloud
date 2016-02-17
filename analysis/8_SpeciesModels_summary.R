# Load the libraries and set working directory
source("analysis/setup.R")

## extract environmental data
finputs=list.files("output/sdm/",pattern="modelinput.*.nc$",full=T,recursive=T)
finputs=finputs[1:2]
env=hSDM.ncExtract(files=finputs,what="envdata",unscale=T)
idata=hSDM.ncExtract(files=finputs,what="spdata")%>%filter(trials>0)
idata$value=ifelse(idata$presences>1,1,0)
idata$variable="observations"

## rescale CLD data to 0-1
env[,grepl("CLD",colnames(env))]=env[,grepl("CLD",colnames(env))]/100


## create long unscaled environmental data for plotting
envl=melt(hSDM.ncExtract(files=finputs,what="envdata",unscale=F),
          id.vars=c("species","x","y","cell"))

## extract model results and evaluation
fresults=list.files("output/sdm/",pattern="modeloutput.*.nc",full=T,recursive=T)
fresults=fresults[1:4]

eval=hSDM.ncExtract(files=fresults,what="eval")
coef=hSDM.ncExtract(files=fresults,what="coef")
ac=hSDM.ncExtract(files=fresults,what="autocor")
pred=hSDM.ncExtract(files=fresults,what="predictions")%>%filter(!is.na(pred))
pred$modelname2=ifelse(pred$modelname=="Cloud","Cloud Model","Precipitation Model")

## add predicted p(presence) to idata
pred=left_join(pred,select(idata,species,cell,presences,trials), by=c("species", "cell"),copy=T)
## add binary presence/absence field
pred$pa=factor(pred$presences>0,labels=c("Undetected","Present"))


## Write out evaluation summary table
print(
  xtable(
    select(eval,c(species,modelname,nPresence,nTrials,auc,cor,DIC,GlobalGelmanMPSRF,MoransI,GearyC)),
    caption="Evaluation of distribution models using interpolated precipitation or cloud product"), 
  include.rownames = FALSE,"html",file="manuscript/modeltable.html",digits=2)


## combine to a single dataframe
pd=rbind.data.frame(melt(env,id.vars=c("species","x","y","cell"),measure.vars=c("PPTJUL","CLDJUL")),
                    select(pred,c(species,x,y,cell,variable=modelname2,value=pred)),
                    select(idata,c(species,x,y,cell,variable,value)))

## identify regions for each species
spregs=list(
  "Protea_cynaroides"=data.frame(xmin=20,xmax=21,ymin=-34.1,ymax=-33.7),
  "Rupicola_peruvianus"=data.frame(xmin=-78,xmax=-73,ymin=4,ymax=5.5),
  "Lepidocolaptes_lacrymiger"=data.frame(xmin=-78,xmax=-73,ymin=4,ymax=5.5))

## add region flags for easier subsetting
pd=mutate(pd,reg=ifelse(species=="Protea_cynaroides"&x>=20&x<=21.&y>=-34.1&y<=-33.7|
                          species%in%c("Rupicola_peruvianus","Lepidocolaptes_lacrymiger")&x>=-78&x<=-73&y>=4&y<=5.5,1,0))
pred=mutate(pred,reg=ifelse(species=="Protea_cynaroides"&x>=20&x<=21.&y>=-34.1&y<=-33.7|
                              species%in%c("Rupicola_peruvianus","Lepidocolaptes_lacrymiger")&x>=-78&x<=-73&y>=4&y<=5.5,1,0))
idata=mutate(idata,reg=ifelse(species=="Protea_cynaroides"&x>=20&x<=21.&y>=-34.1&y<=-33.7|
                                species%in%c("Rupicola_peruvianus","Lepidocolaptes_lacrymiger")&x>=-78&x<=-73&y>=4&y<=5.5,1,0))

ranges=by(pred,list(pred$species),function(x) c(range(x$x),range(x$y)))


## create cropped "fortified" data for ggplot
fcoast1=fortify(crop(coast,extent(ranges[[2]])))
fcoast2=fortify(crop(coast,extent(ranges[[1]])))

## write out tables for PLoS Bio
write.csv(ac,"output/sdm/Fig5_sdm_autocorrelation.csv",row.names=F)
write.csv(filter(pred,trials>0), "output/sdm/Fig5_predictions.csv",row.names=F)

library(grid) # needed for arrow function
library(gridExtra)
library(ggplot2)
library(scales)

## load species icons
library(png)
iprot <- rasterGrob(readPNG("figure/ProteaCynaroides.png"))
ilacr <- rasterGrob(readPNG("figure/lacrymigerLacrymiger.png"))



predscale=scale_fill_gradientn(values=c(0,.5,1),colours=c('white','darkgreen','green'),na.value="transparent")


tsize=20

blanktheme=theme(legend.position="none")+theme(axis.line=element_blank(),
                                               axis.text.x=element_blank(),
                                               axis.text.y=element_blank(),
                                               axis.ticks=element_blank(),
                                               axis.title.x=element_blank(),
                                               axis.title.y=element_blank(),
                                               legend.position="none",
                                               panel.background = element_rect(fill = 'white', colour = 'black'),
                                               panel.border=element_rect(fill="transparent",linetype = "solid", colour = "black"),
                                               panel.grid.major=element_blank(),
                                               panel.grid.minor=element_blank(),
                                               plot.background=element_blank(),
                                               text = element_text(size=tsize))

tsp1="Protea_cynaroides"
csp1="darkred"

tsp2="Lepidocolaptes_lacrymiger"
csp2="darkblue"

csbs="black"#grey(.7)

ps1=
  ggplot(filter(pred,species==tsp1),aes(x=x,y=y,fill=pred)) + 
  geom_raster(colour=pred) +
  facet_grid(~modelname) +
  predscale+
  coord_equal(ratio=1.2)+ theme(legend.position="none",
                                legend.margin=unit(0, "cm"),
                                legend.key.height = unit(0.5, "cm"),
                                plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
                                #                              axis.title.x=element_blank(),
                                #                                axis.title.y=element_blank(),
                                strip.text.x=element_blank(),
                                strip.background = element_blank(),
                                text = element_text(size=tsize))+
  geom_point(data=filter(idata,trials>1&presences==0&species==tsp1),
             aes(x=x,y=y,fill=1),pch=1,col=csbs,cex=.2,lwd=2,alpha=.1)+
  geom_point(data=filter(idata,trials>1&presences==1&species==tsp1),
             aes(x=x,y=y,fill=1),pch=3,col=csp1,cex=.2,lwd=3,alpha=.8)+
  ylab("Latitude")+xlab("Longitude")+
  annotate("rect",
           xmin=spregs[[tsp1]]$xmin,xmax=spregs[[tsp1]]$xmax,ymin=spregs[[tsp1]]$ymin,ymax=spregs[[tsp1]]$ymax,
           fill=NA,col="black")+
  geom_path(data=fcoast1,aes(x=long,y=lat,group = group,fill=1),lwd=.2)


ps1r1=
  ggplot(filter(pred,modelname=="Cloud"&species=="Protea_cynaroides"),aes(x=x,y=y,fill=pred)) + 
  geom_raster(colour=pred) +
  predscale+blanktheme+
  coord_equal(ratio=1.2)+
  geom_point(data=filter(idata,reg==1&trials>1&presences==0&species==tsp1),
             aes(x=x,y=y,fill=1),pch=1,col=csbs,cex=.5,lwd=1,alpha=.2)+
  geom_point(data=filter(idata,reg==1&trials>1&presences==1&species==tsp1),
             aes(x=x,y=y,fill=1),pch=3,col=csp1,cex=1,lwd=1)+
  ylim(unlist(spregs[[tsp1]][,c("ymin","ymax")]))+  xlim(unlist(spregs[[tsp1]][,c("xmin","xmax")]))

ps1r2=
  ggplot(filter(pred,modelname=="Precipitation"&species==tsp1),aes(x=x,y=y,fill=pred)) + 
  geom_raster(colour=pred) +
  predscale+blanktheme+
  coord_equal(ratio=1.2)+
  geom_point(data=filter(idata,reg==1&trials>1&presences==0&species==tsp1),
             aes(x=x,y=y,fill=1),pch=1,col=csbs,cex=.5,lwd=1,alpha=.2)+
  geom_point(data=filter(idata,reg==1&trials>1&presences==1&species==tsp1),
             aes(x=x,y=y,fill=1),pch=3,col=csp1,cex=1,lwd=1)+
  ylim(unlist(spregs[[tsp1]][,c("ymin","ymax")]))+  xlim(unlist(spregs[[tsp1]][,c("xmin","xmax")]))

## second species
ps2=
  ggplot(filter(pred,species==tsp2),aes(x=x,y=y,fill=pred)) + 
  geom_raster(colour=pred) +
  facet_grid(~modelname) +
  predscale+
  coord_equal(ratio=1.3)+ theme(legend.position="right",                             
                                plot.margin=unit(c(0,0,0,0), "cm"),
                                panel.margin=unit(.75, "cm"),
                                legend.margin=unit(0, "cm"),
                                legend.key.width = unit(0.25, "cm"),
                                axis.title.x=element_blank(),
                                #                                axis.title.y=element_blank(),
                                strip.text.x=element_blank(),
                                strip.background = element_blank(),
                                text=element_text(size=tsize))+
  geom_point(data=filter(idata,trials>1&presences==0&species==tsp2),
             aes(x=x,y=y,fill=1),pch=1,col=csbs,cex=.8,lwd=2,alpha=.3)+
  geom_point(data=filter(idata,trials>1&presences==1&species==tsp2),
             aes(x=x,y=y,fill=1),pch=3,col=csp2,cex=2,lwd=3,alpha=1)+
  annotate("rect",
           xmin=spregs[[tsp2]]$xmin,xmax=spregs[[tsp2]]$xmax,ymin=spregs[[tsp2]]$ymin,ymax=spregs[[tsp2]]$ymax,
           fill=NA,col="black")+
  labs(fill = "p(pres)")+
  ylab("Latitude")+
  geom_path(data=fcoast2,aes(x=long,y=lat,group = group,fill=1),lwd=.5)



ps2r1=
  ggplot(filter(pred,modelname=="Cloud"&species==tsp2),aes(x=x,y=y,fill=pred)) + 
  geom_raster(colour=pred) +
  predscale+blanktheme+
  coord_equal(ratio=1.3)+
  geom_point(data=filter(idata,reg==1&trials>1&presences==0&species==tsp2),
             aes(x=x,y=y,fill=1),pch=1,col=csbs,cex=1,lwd=2,alpha=.3)+
  geom_point(data=filter(idata,reg==1&trials>1&presences==1&species==tsp2),
             aes(x=x,y=y,fill=1),pch=3,col=csp2,cex=2,lwd=3)+
  ylim(unlist(spregs[[tsp2]][,c("ymin","ymax")]))+  xlim(unlist(spregs[[tsp2]][,c("xmin","xmax")]))+
  geom_path(data=fcoast2,aes(x=long,y=lat,group = group,fill=1),lwd=.2)

ps2r2=
  ggplot(filter(pred,modelname=="Precipitation"&species==tsp2),aes(x=x,y=y,fill=pred)) + 
  geom_raster(colour=pred) +
  predscale+blanktheme+
  coord_equal(ratio=1.3)+
  geom_point(data=filter(idata,reg==1&trials>1&presences==0&species==tsp2),
             aes(x=x,y=y,fill=1),pch=1,col=csbs,cex=1,lwd=2,alpha=.3)+
  geom_point(data=filter(idata,reg==1&trials>1&presences==1&species==tsp2),
             aes(x=x,y=y,fill=1),pch=3,col=csp2,cex=2,lwd=3)+
  ylim(unlist(spregs[[tsp2]][,c("ymin","ymax")]))+  xlim(unlist(spregs[[tsp2]][,c("xmin","xmax")]))+
  geom_path(data=fcoast2,aes(x=long,y=lat,group = group,fill=1),lwd=.2)



## autocorrelation
pac1=
  ggplot(ac, aes(x=dist2, y=mean,group=interaction(modelname,species),linetype=modelname,colour=species))+
  geom_linerange(aes(ymax = mean+sd, ymin=mean-sd),linetype="solid",alpha=.2,shape=NA)+
  #geom_point(alpha=.5)+
  geom_line(lwd=.5)+
  geom_hline(yintercept=0)+
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  ylab("Spatial\nAutocorrelation")+xlab("Distance (km)")+
  theme(legend.position="none",
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0,0,0), "cm"),
        text=element_text(size=tsize))+
  annotation_logticks(sides = "b")+
  scale_y_continuous(breaks=seq(0, 1, .5))+
  scale_colour_manual(name = "Species",
                      values = c(csp2,csp1)) +   
  scale_shape_manual(name = "Species",
                     values = c(1,3)) +   
  scale_linetype_manual(name = "Model",
                        values = c("solid","dashed"))+
  guides(colour = guide_legend(nrow=1),
         linetype = guide_legend(nrow=1))

## boxplots
pbox=
  ggplot(filter(pred,pred$trials>=5), aes(x=pa=="Present", y=pred,fill=species,linetype=modelname))+ 
  theme(legend.position="none",
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0,0,0), "cm"),
        #        panel.margin=unit(c(0,0,0,0), "cm"),
        #        axis.title.x=element_blank(),
        text=element_text(size=tsize),
        #        axis.text.x = element_text(angle = 20, hjust = 1),
        strip.background = element_blank(), strip.text = element_blank())+
  scale_fill_manual(name = "Model",values = c(csp2,csp1)) +   
  geom_boxplot(notch=T,varwidth=T,outlier.colour=grey(.8),outlier.size=.8,position=position_dodge(.8),lwd=.8,col="black")+
  facet_wrap(~species)+
  scale_y_continuous(breaks=seq(0, 1, .5))+
  xlab("Species Observed")+ylab("p(presence)")


## make single plot with two species and autocorrelation
pdf(file=paste0("figure/SDM_overview.pdf"),
    width=10,height=16,pointsize=24)

png(file=paste0("figure/SDM_overview.png"),
    width=2800,height=4100,pointsize=24,res=300)

print(ps1, vp = viewport(width = 1, height = .3, x=.5,y=.375))
print(ps1r1, vp = viewport(width = .35, height = .25, x = .385, y = 0.42))
print(ps1r2, vp = viewport(width = .35, height = .25, x = .85, y = 0.42))

print(ps2, vp = viewport(width = 1, height = .55, x=.55,y=.74))
print(ps2r1, vp = viewport(width = .37, height = .25, x = .38, y = 0.8))
print(ps2r2, vp = viewport(width = .37, height = .25, x = .75, y = 0.8))

## topleft
pushViewport(viewport(width = .09, height = .09, x=.21,y=.97))
print(grid.draw(ilacr))
popViewport()
## top right
pushViewport(viewport(width = .09, height = .09, x=.57,y=.97))
print(grid.draw(ilacr))
popViewport()
## top left
pushViewport(viewport(width = .07, height = .07, x=.14,y=.44))
print(grid.draw(iprot))
popViewport()
# top right
pushViewport(viewport(width = .07, height = .07, x=.6,y=.44))
print(grid.draw(iprot))
popViewport()

### bottom row of plots (ac and boxplot)
print(pac1, vp = viewport(width = .45, height = .24, x = .25, y = 0.122))
print(pbox, vp = viewport(width = .5, height = .25, x = .75, y = 0.129))

## ac plot
pushViewport(viewport(width = .09, height = .09, x=.42,y=.21))
print(grid.draw(ilacr))
popViewport()

pushViewport(viewport(width = .09, height = .09, x=.22,y=.15))
print(grid.draw(iprot))
popViewport()

## boxplot

pushViewport(viewport(width = .09, height = .09, x=.65,y=.21))
print(grid.draw(ilacr))
popViewport()

pushViewport(viewport(width = .09, height = .09, x=.85,y=.21))
print(grid.draw(iprot))
popViewport()


## add panel labels
pushViewport(viewport())
tgp=gpar(cex = .5, col = "black")
grid.text(label = "a" ,x = 0.06,y = .98,gp = tgp)
grid.text(label = "b" ,x = 0.52,y = .98,gp = tgp)
grid.text(label = "c" ,x = 0.06,y = .48,gp = tgp)
grid.text(label = "d" ,x = 0.55,y = .48,gp = tgp)
grid.text(label = "e" ,x = 0.06,y = .24,gp = tgp)
grid.text(label = "f" ,x = 0.55,y = .24,gp = tgp)
popViewport()

dev.off()


## Model input data
penv1=ggplot(filter(envl,species=="Protea_cynaroides"),aes(x=x,y=y,fill=value)) + 
  geom_raster() +
  facet_wrap(~variable,nrow=2) +
  scale_fill_gradientn(colours=c('blue','white','red'),breaks=c(-3,0,3,6),
                       na.value="transparent")+
  coord_equal(ratio=1.2)+ 
  theme(legend.position="bottom",
        legend.margin=unit(0, "cm"),
        legend.key.height = unit(0.5, "cm"))+
  ylab("")+xlab("")+labs(fill = "Standardized Value")+
  geom_path(data=fcoast1,aes(x=long,y=lat,group = group,fill=1),lwd=.2)

penv2=ggplot(filter(envl,species=="Lepidocolaptes_lacrymiger"),aes(x=x,y=y,fill=value)) + 
  geom_raster() +
  facet_wrap(~variable,nrow=2) +
  scale_fill_gradientn(colours=c('blue','white','red'),breaks=c(-3,0,3,6),
                       na.value="transparent",name="Standardized Value")+
  coord_equal(ratio=1.2)+ 
  theme(legend.position="none",
        legend.margin=unit(0, "cm"),
        legend.key.height = unit(0.5, "cm"))+
  ylab("")+xlab("")+labs(fill = "Standardized Value")+
  geom_path(data=fcoast2,aes(x=long,y=lat,group = group,fill=1),lwd=.2)

png(file=paste0("figure/SDM_inputs.png"),width=3000,height=3000,pointsize=38,res=300)
print(penv1,vp = viewport(width = 1, height = .5, x=.5,y=.25))
print(penv2,vp = viewport(width = 1, height = .6, x=.5,y=.75))
dev.off()

#### Separate figures for poster
png(file=paste0("figure/SDM_overview1_CFR.png"),  width=4000,height=4000,pointsize=38,res=1200)
print(ps1)
print(ps1r1, vp = viewport(width = .35, height = .25, x = .385, y = 0.42))
print(ps1r2, vp = viewport(width = .35, height = .25, x = .85, y = 0.42))
dev.off()

png(file=paste0("figure/SDM_overview1_SouthAmerica.png"),
    width=4000,height=4000,background="transparent",pointsize=38,res=1200)
print(ps2)
dev.off()

png(file=paste0("figure/SDM_overview1_CFR.png"),width=3000,height=3000,pointsize=38,res=600)
###  (ac and boxplot)
print(pac1)
print(pbox)
dev.off()


## Range size statistics
rs=pred%.%
  group_by(species,modelname)%.%
  summarise(area=sum(pred)) %>%
  spread(key="modelname", value=area) %>%
  mutate(pchange=(Precipitation-Cloud)/Precipitation)

rs



######################################
### Old stuff below here

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

png(file=paste0("figure/SDM_all.png"),width=3000,height=3000,pointsize=48)
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


## predictions with expert range
pred$type=paste(pred$species,pred$modelname2)

ggplot(filter(pred),aes(x=x,y=y,fill=pred)) + geom_tile() +
  facet_wrap(~type,scales="free") +
  scale_fill_gradientn(values=c(0,.5,1),colours=c('white','blue','red'),na.value="transparent")+
  coord_equal(ratio=1)+ theme(legend.position="right")+
  geom_point(data=filter(idata,presences==1&species=="Lepidocolaptes_lacrymiger"),aes(fill=1),pch=3,cex=3,lwd=3)+
  ylab("")+xlab("")+
  labs(fill = "")+
  theme(text=element_text(size=12))


## compare predictions
p1=
  gplot(filter(pred,aes(x=x,y=y,fill=pred,maxpixels=1e4) + 
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
                 #geom_pointrange(aes(ymax = max, ymin=min, group=modelname,colour=modelname))+
                 geom_pointrange(aes(ymax = mean+sd, ymin=mean-sd, group=modelname,colour=modelname))+
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
               
               