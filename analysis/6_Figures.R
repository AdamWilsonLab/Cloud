### Figures and tables for MOD09 Cloud Manuscript
source("analysis/setup.R")
source("analysis/6_LoadData.R")

n=100
res=1e7
greg=list(ylim=c(-60,84),xlim=c(-180,180))
    
## Figure 1: 4-panel summaries
#- Annual average
n=100
colR(n)
cwidth=.5
#quantile(cf_mean,c(0,0.01,0.1,0.25,.50,.75,.90,.99,1),na.rm=T)

p_mean=  levelplot(cf_mean,col.regions=bgr(1:100,n,br=50.9)$col,cuts=99,at=seq(0,100,len=n),
                 colorkey=list(space="left",width=cwidth,height=.5,labels=list(at=c(0, 50, 100)),border=NA),
                 scales=list(draw=F),
                 panel=panel.levelplot.raster,margin=F,maxpixels=res,ylab="",xlab="",useRaster=T,ylim=greg$ylim)+
    latticeExtra::layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col=grey(.2)),under=T)+
  latticeExtra::layer(sp.lines(coast,lwd=.5,),under=F)
 p_mean$strip=strip.custom(factor.levels="a. Mean Cloud Frequency (%)")
 p_mean$par.strip.text=list(cex=1)
#quantile(spatial,c(0,0.01,0.1,0.25,.50,.75,.90,.99,1),na.rm=T)
#0%    1%   10%   25%   50%   75%   90%   99%  100% 
#0.01  0.38  0.50  0.65  1.05  2.42  4.99 11.69 57.50 

#c_spatial=bgr(log10(seq(.01,7,len=n)),n,br=log10(0.45))
c_spatial=bgr(seq(0,25,len=n),n,br=3,c1=c("darkblue","blue","grey"))

p_spatial=levelplot(spatial,margin=F,maxpixels=res,
                    col.regions=c_spatial$col,at=c_spatial$at,
                    panel=panel.levelplot.raster,ylim=greg$ylim,ylab="",xlab="",#zscaleLog=T,
                    colorkey=list(title="SD Cloud Frequency", space="left",width=cwidth,height=.5,
                                  labels=list(
                    at=c(0,10,20),
                    labels= c(0,10,20))),
scales=list(draw=F),useRaster=T)+
  latticeExtra::layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col=grey(.2)),under=T)+
  latticeExtra::layer(sp.lines(coast,col="black",lwd=.5),under=F)
p_spatial$strip=strip.custom(factor.levels="c. Spatial Variability (SD)") #  latticeExtra::layer(sp.lines(coast,lwd=.5,),under=F)
p_spatial$par.strip.text=list(cex=1)
                    
# quantile(intra,c(.50,.75,.90,.99,1),na.rm=T)
#50%   75%   90%   99%  100% 
#10.59 15.26 23.11 36.30 98.43 


c_intra=bgr(log10(seq(1,37,len=n)),n,br=log10(10.59))
c_intra=bgr(seq(0,40,len=n),n,br=10.59)

#c_intra=bgr(seq(0,25,len=n),n,br=10.59)

# p_intra=levelplot(intra,
#                   col.regions=c_intra$col,at=c_intra$at,cuts=n-1,
#                   margin=F,maxpixels=res,ylab="",xlab="",
#                   panel=panel.levelplot.raster,ylim=greg$ylim,
#                   colorkey=list(space="right",height=.5,width=cwidth,
# #                                labels=list(at=log10(c(1,5,10,36)),
# #                                            labels=c(1,5,10,36))),
#                                 labels=list(at=c(0,20,40),labels=c(0,20,40))),                  
#                   useRaster=T,scales=list(draw=F))+
#   latticeExtra::layer(panel.polygon(x=c(-180,-180,180,180),y=c(-61,86,86,-61),col=grey(.2)),under=T)+
#   latticeExtra::layer(sp.lines(coast,col="black",lwd=.5),under=F)
# p_intra$strip=strip.custom(factor.levels="b. Intra-annual Variability (SD)") 
# p_intra$par.strip.text=list(cex=1)
#                   
# #quantile(inter,c(.50,.75,.90,.99,1),na.rm=T)
# #50%    75%    90%    99%   100% 
# #11.39  13.31  16.20  20.61 140.77 
# 
# #c_inter=bgr(log10(seq(1,21,len=n)),n,br=log10(11.39))
# c_inter=bgr(seq(0,25,len=n),n,br=11.39)
# 
# # 
# p_inter=levelplot(inter,cuts=n-1,margin=F,maxpixels=res,
#                   col.regions=c_intra$col,at=c_intra$at,
#                   panel=panel.levelplot.raster,ylim=greg$ylim,ylab="",xlab="",#zscaleLog=10,
#                   colorkey=list(title="Cloud Frequency (%)", space="right",height=.5,width=cwidth,
# #                                labels=list(at=log10(c(1,seq(2,8,by=2),10,seq(20,80,by=20),100)),
# #                                            labels=c(1,rep("",4), 10,rep("",4), 100))),useRaster=T,
# #                  labels=list(at=log10(c(1,5,10,21)),
# #                              labels=c(1,5,10,21))),useRaster=T,
#                   labels=list(at=c(0,10,20)),
#                               labels=c(1,10,20)),useRaster=T,
# 
#                   scales=list(draw=F))+
# #                  main=textGrob("    d       Inter-annual Variability (SD)", x = 0, hjust = 0,gp=gpar(fontface="bold")))+
#   latticeExtra::layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col=grey(.2)),under=T)+
#   latticeExtra::layer(sp.lines(coast,col="black",lwd=.5),under=F)
# p_inter$strip=strip.custom(factor.levels="d. Inter-annual Variability (SD)") #  latticeExtra::layer(sp.lines(coast,lwd=.5,),under=F)
# p_inter$par.strip.text=list(cex=1)
svar=stack(inter,intra)
NAvalue(svar)=255
# 
p_var=levelplot(svar,cuts=n-1,margin=F,maxpixels=res,
                  col.regions=c_intra$col,at=c_intra$at,
                  panel=panel.levelplot.raster,ylim=greg$ylim,ylab="",xlab="",#zscaleLog=10,
                  colorkey=list(title="Cloud Frequency (%)", space="right",
                                height=.5,width=cwidth,
                                labels=list(at=c(0,20,40)),
                                labels=c(1,20,40)),useRaster=T,
                  
                  scales=list(draw=F),layout=c(1,2))+
  latticeExtra::layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col=grey(.2)),under=T)+
  latticeExtra::layer(sp.lines(coast,col="black",lwd=.5),under=F)
p_var$strip=strip.custom(factor.levels=c("b. Inter-annual Variability (SD)","d. Intra-annual Variability (SD)")) #  latticeExtra::layer(sp.lines(coast,lwd=.5,),under=F)


ndifx=.06
ndify=.06

#pdf("manuscript/figures/MeanInter2.pdf",width=8,height=4,pointsize=24,bg="transparent")
png("manuscript/figures/MeanInter2.png",width=8,height=4,units="in",res=1200,pointsize=24,bg="white")
#trellis.par.set(my.theme)
print(p_mean,position=c(0,.5-ndify,.5+ndifx,1),more=T)
#print(p_intra,position=c(.5-ndify,.5-ndifx,1,1),more=T)
print(p_spatial,position=c(0,0,.5+ndifx,.5+ndify),more=T)
#print(p_inter,position=c(.5-ndifx,0,1,.5+ndify),more=F)
print(p_var,position=c(.5-ndifx,0,1,1),more=F)
                         
dev.off()

## Scatterplots of variables
png("manuscript/figures/scatterplot.png",width=2100,height=2000,res=600,pointsize=42,bg="white")
trellis.par.set(my.theme)
splom(stack(cf_mean,inter,intra,spatial),maxpixels=1e6)
dev.off()


#########
#### Seasonality plot

col=read.csv("data/MCD09_deriv/coltable.csv",stringsAsFactors=F)



pCirc=function(r,n=100) {
  angs=seq(0,360,len=100)
  #lapply(r,function(r) {
  cbind(
    y=r*cos((pi/180)*(-angs)),
    x=r*sin((pi/180)*(-angs))
  )
  #})
}


mangle=(-seq(30,360,30)+prot+15)*(pi/180)
rmons=c(1:3,9:12)
lmons=4:8
mon=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec") # month names
ladj=65
lims=c(-80,80)


####################################
## Seasonality plot

## color key
k1=xyplot(y~x,col=col$val[col$exists],data=col[col$exists,],pch=16,cex=1.2,
       xlab="",ylab="",scales=list(draw=F),ylim=lims,xlim=lims,asp=1,
#          main="Seasonal Cloud Concentration"
       par.settings = list(axis.line = list(col = "transparent")))+
  latticeExtra::layer(panel.polygon(pCirc(r=20,n=100),border="grey"))+
  latticeExtra::layer(panel.polygon(pCirc(r=40,n=100),border="grey"))+
  latticeExtra::layer(panel.polygon(pCirc(r=60,n=100),border="grey"))+
  latticeExtra::layer(panel.segments(0,0,60*cos(mangle-(15*pi/180)),60*sin(mangle-(15*pi/180)),col="grey"))+ #draw angles
  #latticeExtra::layer(panel.text(x=ladj*cos(mangle[lmons]),y=ladj*sin(mangle[lmons]),mon[lmons],
  #                               pos=4,cex=1,offset=0,srt=(mangle[lmons])*180/pi))+ #add left months))
  #latticeExtra::layer(panel.text(x=ladj*cos(mangle[rmons]),y=ladj*sin(mangle[rmons]),mon[rmons],
  #                               pos=2,cex=1,offset=0,srt=((mangle[rmons])*180/pi)-180))+ # add right months
  latticeExtra::layer(panel.text(x=ladj*cos(mangle[lmons[-1]]),y=ladj*sin(mangle[lmons[-1]]),mon[lmons[-1]],
                                 pos=4,cex=1,offset=0,srt=(mangle[lmons[-1]])*180/pi))+ #add left months))
  latticeExtra::layer(panel.text(x=ladj*cos(mangle[rmons[c(1,2,6,7)]]),y=ladj*sin(mangle[rmons[c(1,2,6,7)]]),mon[rmons[c(1,2,6,7)]],
                                 pos=2,cex=1,offset=0,srt=((mangle[rmons[c(1,2,6,7)]])*180/pi)-180))+ # add right months
    latticeExtra::layer(panel.text(x=-2,y=c(-20,-40,-60),c(20,40,60),pos=4,cex=1,col=c("white","black","black"))) #add scale text


g1=levelplot(cf_visseas,col.regions=cf_visseas@legend@colortable,cuts=length(cf_visseas@legend@colortable),at=0:length(cf_visseas@legend@colortable),
             colorkey=F,panel=panel.levelplot.raster,margin=F,maxpixels=res,ylab="",xlab="",useRaster=T,ylim=c(-60,70),
             asp=1,scales=list(draw=F))+#cex=1,y=list(at=c(-40,0,40))))+
  latticeExtra::layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col=grey(.2)),under=T)+
  latticeExtra::layer(sp.polygons(as(regs[["CFR2"]],'SpatialPolygons'),col="red",lwd=1.5),under=F)+
  latticeExtra::layer(sp.polygons(as(regs[["Venezuela2"]],'SpatialPolygons'),col="red",lwd=1.5),under=F)+
  latticeExtra::layer(sp.lines(coast,col="white",lwd=.5),under=F)

## regional plots
r_cfr=crop(cf_visseas,regs[["CFR2"]],datatype="INT2U")
r1=levelplot(r_cfr,col.regions=cf_visseas@legend@colortable,
             cuts=length(cf_visseas@legend@colortable),
             at=0:length(cf_visseas@legend@colortable),
             colorkey=F,panel=panel.levelplot.raster,margin=F,
             maxpixels=res,ylab="",xlab="",useRaster=T,asp=1,
             scales=list(cex=1,y=list(at=c(-34,-30,-26))))+
  latticeExtra::layer(sp.lines(coast,col="white",lwd=.5),under=F)

r_ven=crop(cf_visseas,regs[["Venezuela2"]],datatype="INT2U")
r2=levelplot(r_ven,col.regions=cf_visseas@legend@colortable,
             cuts=length(cf_visseas@legend@colortable),
             at=0:length(cf_visseas@legend@colortable),
             asp=1,colorkey=F,panel=panel.levelplot.raster,margin=F,
             maxpixels=res,ylab="",xlab="",useRaster=T,asp=1,
             scales=list(cex=1,y=list(at=c(2,6,10))))+
  latticeExtra::layer(sp.lines(coast,col="white",lwd=.5),under=F)

## draw it
png(width=2400,height=1600,res=300,pointsize=16,
    type="cairo-png",file="manuscript/figures/Seasonality.png")
print(g1,position=c(0.4,.5,1,1),more=T) #global
print(r1,position=c(.51,-0.05,1,.7),more=T)
print(r2,position=c(0,-.05,.55,.7),more=T)
print(k1,position=c(-.10,.4,.55,1.1),just=c("left","top"),more=T) #legend
## add panel labels
gp=gpar(fontsize=16, col="black",fontface="bold")
grid.text("a",x=0.01,y=.995,just=c("left","top"),gp=gp)
grid.text("b",x=0.42,y=.995,just=c("left","top"),gp=gp)
grid.text("c",x=0.05,y=.52,just=c("left","top"),gp=gp)
grid.text("d",x=0.57,y=.52,just=c("left","top"),gp=gp)
dev.off()


## Draw separate figures for key and full map
pdf("manuscript/figures/SeasonalityMap.pdf",width=16,height=8,pointsize=24,bg="transparent")
print(g1) #global
dev.off()

pdf("manuscript/figures/SeasonalityKey.pdf",width=4,height=4,pointsize=48,bg="transparent")
#png("manuscript/figures/SeasonalityKey.png",width=2000,height=2000,res=1200,pointsize=48,bg="transparent")
print(k1) #legend
dev.off()



####################################################################
### Regional Comparisons
## Compare with worldclim and NPP
#wc=stack(as.list(paste("/mnt/data/jetzlab/Data/environ/global/worldclim/prec_",1:12,".bil",sep="")))
wc_map=stack(as.list(paste("/mnt/data/jetzlab/Data/environ/global/worldclim/bio_12.bil",sep="")))
wc_dem=stack(as.list(paste("/mnt/data/jetzlab/Data/environ/global/worldclim/alt.bil",sep="")))


## read in GEWEX 1-degree data
gewex=mean(brick("data/src/gewex/CA_PATMOSX_NOAA.nc",varname="a_CA"))
names(gewex)="PATMOS-x GEWEX AVHRR"

## calculate 1-degree means of MODCF data
#MOD_gewex=gewex
#MOD_gewex@data@values=1:length(MOD_gewex@data@values)
#MOD_gewex2=zonal(mod09a,MOD_gewex,fun='mean')
res=1e4

r="Venezuela2"

cf_r=crop(cf_mean,regs[[r]])
tmap=crop(wc_map,regs[[r]])


# ylab.right = "Cloud Frequency (%)",par.settings = list(layout.widths = list(axis.key.padding = 0.1,axis.left=0.6,ylab.right = 3,right.padding=2)),
pars=list(layout.heights=list(key.bottom=2,key.top=1),layout.widths = list(axis.key.padding = 3,axis.left=1))
p1=levelplot(cf_r,col.regions=colR(n),at=seq(0,100,len=99),
             ylab="",xlab="", scales=list(draw=F),
             colorkey=list(space="left",width=1,height=.75,labels=list(labels=c(0,50,100),at=c(0,50,100))),
    cuts=99,margin=F,maxpixels=1e7,par.settings = pars)+
  latticeExtra::layer(sp.lines(coast,col="black",lwd=.5),under=F)

p2=levelplot(crop(gewex,regs[[r]]),col.regions=colR(n),at=seq(0,1,len=99),cuts=99,margin=F,max.pixels=1e6,
             ylab="",xlab="",
             colorkey=F,scales=list(y=list(draw=F),x=list(draw=T)),
#             colorkey=list(space="right",width=1,height=.75,labels=list(labels=c(0,50,100),at=c(.0,.5,1))),
    par.settings = pars)+
  latticeExtra::layer(sp.lines(coast,col="black",lwd=.5),under=F)

p3=levelplot(tmap,col.regions=rev(terrain.colors(n)),cuts=100,
             at=seq(tmap@data@min,tmap@data@max,len=100),margin=F,maxpixels=1e6,
             scales=list(draw=F),
             colorkey=list(space="right",height=.5,width=1,
                           labels=list(labels=c(2000,6000,10000),at=c(2000,6000,10000))),
                           xlab="",ylab="",
             useRaster=T,
    par.settings = pars)+
  latticeExtra::layer(sp.lines(coast,col="black",lwd=.5),under=F)

p4=levelplot(crop(wc_dem,regs[[r]]),col.regions=terrain.colors(n),
             cuts=99,margin=F,max.pixels=1e6,
             ylab="",xlab="",maxpixels=1e7,
#             scales=list(draw=F),
             colorkey=list(space="right",height=.5,width=1,
                           labels=list(labels=c(0,2500,5000),
                           at=c(0,2500,5000))),
    par.settings = pars,scales=list(y=list(at=c(2,6,10))))+
  latticeExtra::layer(sp.lines(coast,col="black",lwd=.5),under=F)

## combine
lside=c("a    MODCF (%)"=p1,"c     PATMOS-x GEWEX (%)"=p2,
        layout=c(1,2),merge.legends=T,x.same=T,y.same=T)
rside=c("b    WorldClim Precip (mm)"=p3,"d    Elevation (m)"=p4,
        layout=c(1,2),x.same=T,y.same=T,merge.legends=T)
      

png("manuscript/figures/Resolution2.png",width=3000,height=1500,res=300,pointsize=47,bg="white")
trellis.par.set(my.theme)
#pdf("output/mod09_resolution.pdf",width=11,height=8.5)
#print(c("b    WorldClim Precip (mm)"=p3,"d    Elevation (m)"=p4,x.same=T,y.same=T,merge.legends=T,layout=c(1,2))
      #print(c("a     PATMOS-x GEWEX (%)"=p2,"b    MODCF (%)"=p1,"c    Elevation (m)"=p4,x.same=T,y.same=T,merge.legends=T,layout=c(1,3)))
print(lside,position=c(0,0,.5,1),more=T)
print(rside,position=c(0.48,0,1,1),more=F)
dev.off()

###  Spatial Autocorrelation

l1=xyplot(conc~theta,col=col$val[col$exists],data=col[col$exists,],pch=16,cex=1,
          scales=list(
            x=list(labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
                   at=seq(15,360,30))),
          xlab="Month",ylab="Concentration (%)")
l1b=xyplot(conc~theta,col=col$rgbcol,data=col,pch=16,cex=1,
          scales=list(
            x=list(labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
                   at=seq(15,360,30))),
          xlab="Month",ylab="Concentration (%)")
c(l1,l1b)

l2=levelplot(x~x*y,colours=seasl$col,data=seasl,pch=16,cex=1,
            panel = function(x, y, colours, subscripts, ...) {
            panel.xyplot(x, y, pch = 21, col = "transparent",fill = colours[subscripts])
            },ylab="",xlab="",
          colorkey=F,scales=list(draw=T))+
    layer(sp.lines(coast,col="black"),under=F)

pdf(width=11,height=8.5,file="manuscript/figures/Seasonality.pdf",useDingbats=F,pointsize=12)
trellis.par.set(my.theme)

print(l1,position=c(0,0,1,.25),more=T)
print(l2,position=c(0,.2,1,1),more=F)

dev.off()



##############
### validation stations
v1=xyplot(lat~lon,groups=era, data=st@data,asp=.6,auto.key=F,
          ylab="Latitude",xlab="Longitude",ylim=greg$ylim,
          par.settings=list(superpose.symbol =list(cex=.5,pch=16,col=c("black","red"))))+
  latticeExtra::layer(sp.lines(coast,col="black",lwd=1),under=T)

png("manuscript/figures/ValidationStations.png",width=2000,height=1500,res=300,pointsize=47,bg="white")
trellis.par.set(my.theme)
print(v1)
dev.off()

##########################
## Sahara orbital correction
dt=stack("data/out/SaharaBiasCorrectionExample.tif")[[1:2]]
names(dt)=c("a   Uncorrected_Terra ","b    Corrected_Terra")
gain(dt)=.01

pars=list(layout.heights=list(key.bottom=2,key.top=1),layout.widths = list(axis.key.padding = 3,axis.left=1))
bcols=colorRampPalette(c("#08306b","#0d57a1","#2878b8","#4997c9","#72b2d7","#a2cbe2","#c7dcef","#deebf7","#f7fbff"))

png("manuscript/figures/Sahara.png",width=2000,height=1500,res=300,pointsize=47,bg="white")
trellis.par.set(my.theme)
s1=levelplot(dt[[1]],col.regions=bcols(100),at=seq(0,75,len=101),
          cuts=99,margin=F,max.pixels=1e6,par.settings = pars)+layer(sp.lines(coast))
s2=levelplot(dt[[2]],col.regions=bcols(100),at=seq(0,75,len=101),
          cuts=99,margin=F,max.pixels=1e6,par.settings = pars)+layer(sp.lines(coast))
print(c("a   Uncorrected Terra"=s1,"b    Corrected Terra"=s2))
dev.off()

## vector plot of seasonality
vectorplot(seas,narrows=2e3, lwd.arrows=0.6, length=unit(5e-2, 'npc'),
           maxpixels=1e5, region=TRUE, margin=FALSE,
           isField=TRUE, reverse=FALSE,
           unit='degrees', scaleSlope=TRUE,
           aspX=0.08)
