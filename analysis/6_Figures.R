### Figures and tables for MOD09 Cloud Manuscript
source("analysis/setup.R")
source("analysis/6_LoadData.R")

n=100
res=1e6
greg=list(ylim=c(-60,84),xlim=c(-180,180))
    
## Figure 1: 2-panel summaries
#- Annual average
n=100
colR(n)
#cellStats(cf_mean,median)
p_mean=levelplot(cf_mean,col.regions=bgr(1:100,n,br=50.9)$col,cuts=99,at=seq(0,100,len=n),
                 colorkey=list(space="right",height=.75),
                 main=textGrob("        a       Mean Cloud Frequency (%)", x = 0, hjust = 0,gp=gpar(fontface="bold")),
                 scales=list(x=list(draw=F)),
          panel=panel.levelplot.raster,margin=F,maxpixels=res,ylab="",xlab="",useRaster=T,ylim=greg$ylim)+
#    layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col="white"),under=T)#+
    layer(sp.lines(coast,lwd=.5,),under=F)

#p_intra=levelplot(intra,col.regions=bgr(1:55,n,br=8)$col,cuts=99,at=seq(0,55,length=n),margin=F,maxpixels=res,
#                  panel=panel.levelplot.raster,ylim=greg$ylim,colorkey=list(space="right",height=.75),useRaster=T)+
#    layer(panel.polygon(x=c(-180,-180,180,180),y=c(-61,86,86,-61),col="black"),under=T)#+
#    layer(sp.lines(coast,col="black",lwd=.5),under=F)

#
# cellStats(inter,quantile,c(.50,.75,.90,.99,1),na.rm=T)

p_inter=levelplot(inter,col.regions=bgr(1:30,n,br=11)$col,cuts=99,at=seq(0,30,length=100),margin=F,maxpixels=res,
                  panel=panel.levelplot.raster,ylim=greg$ylim,ylab="",xlab="",
                  colorkey=list(title="Cloud Frequency (%)", space="right",height=.75),useRaster=T,
                  main=textGrob("        b       Inter-annual Variability (SD)", x = 0, hjust = 0,gp=gpar(fontface="bold")))+
  #    layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col="black"),under=T)#+
    layer(sp.lines(coast,col="black",lwd=.5),under=F)


#pdf("manuscript/figures/Figures.pdf",width=11,height=8.5,pointsize=14)
png("manuscript/figures/MeanInter.png",width=2100,height=2000,res=300,pointsize=42,bg="white")
trellis.par.set(my.theme)
#p3=c("Mean Cloud Frequency (%)"=p_mac,"Max Cloud Frequency (%)"=p_max,"Interannual Variability (sd)"=p_inter,"Intraannual Variability (sd)"=p_intra,x.same=T,y.same=F,merge.legends=T,layout=c(2,2))
#p3=c("Mean Cloud Frequency (%)"=p_mean,"Inter-annual Variability (SD)"=p_inter,"Intra-annual Variability (SD)"=p_intra,x.same=T,y.same=T,merge.legends=T,layout=c(1,3))
print(p_mean,position=c(0,.53,1,1),more=T)
print(p_inter,position=c(0,0,1,.52),more=F)
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

png(width=1000,height=1000,pointsize=38,file="manuscript/figures/SeasKey_%0d.png")
### Create the color key
xyplot(conc~theta,col=col$val[col$exists],data=col[col$exists,],pch=16,cex=1,
       scales=list(
         x=list(labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
                at=seq(15,360,30))),
       xlab="Month",ylab="Concentration (%)")
dev.off()

####################################
## Seasonality plot

## color key
k1=xyplot(y~x,col=col$val[col$exists],data=col[col$exists,],pch=16,cex=1.2,
       xlab="",ylab="",scales=list(draw=F),ylim=lims,xlim=lims,asp=1,
#          main="Seasonal Cloud Concentration"
       par.settings = list(axis.line = list(col = "transparent")))+
  layer(panel.polygon(pCirc(r=20,n=100),border="grey"))+
  layer(panel.polygon(pCirc(r=40,n=100),border="grey"))+
  layer(panel.polygon(pCirc(r=60,n=100),border="grey"))+
  layer(panel.segments(0,0,60*cos(mangle-(15*pi/180)),60*sin(mangle-(15*pi/180)),col="grey"))+ #draw angles
  layer(panel.text(x=ladj*cos(mangle[lmons]),y=ladj*sin(mangle[lmons]),mon[lmons],pos=4,cex=1,offset=0,srt=(mangle[lmons])*180/pi))+ #add left months))
  layer(panel.text(x=ladj*cos(mangle[rmons]),y=ladj*sin(mangle[rmons]),mon[rmons],pos=2,cex=1,offset=0,srt=((mangle[rmons])*180/pi)-180))+ # add right months
  layer(panel.text(x=-2,y=c(-20,-40,-60),c(20,40,60),pos=4,cex=1,col=c("white","black","black"))) #add scale text


g1=levelplot(cf_visseas,col.regions=cf_visseas@legend@colortable,cuts=length(cf_visseas@legend@colortable),at=0:length(cf_visseas@legend@colortable),
             colorkey=F,panel=panel.levelplot.raster,margin=F,maxpixels=res,ylab="",xlab="",useRaster=T,ylim=c(-60,70),
             scales=list(cex=1,y=list(at=c(-40,0,40))))+
  layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col="black"),under=T)+
  layer(sp.polygons(as(regs[["CFR2"]],'SpatialPolygons'),col="red",lwd=1.5),under=F)+
  layer(sp.polygons(as(regs[["Venezuela2"]],'SpatialPolygons'),col="red",lwd=1.5),under=F)+
  layer(sp.lines(coast,col="white",lwd=.5),under=F)

## regional plots
r_cfr=crop(cf_visseas,regs[["CFR"]],datatype="INT2U")
r1=levelplot(r_cfr,col.regions=cf_visseas@legend@colortable,cuts=length(cf_visseas@legend@colortable),at=0:length(cf_visseas@legend@colortable),
                colorkey=F,panel=panel.levelplot.raster,margin=F,maxpixels=res,ylab="",xlab="",useRaster=T,scales=list(cex=1,y=list(at=c(-34,-30,-26))),asp=1)+
  layer(sp.lines(coast,col="white",lwd=.5),under=F)
r_ven=crop(cf_visseas,regs[["Venezuela2"]],datatype="INT2U")
r2=levelplot(r_ven,col.regions=cf_visseas@legend@colortable,cuts=length(cf_visseas@legend@colortable),at=0:length(cf_visseas@legend@colortable),
                colorkey=F,panel=panel.levelplot.raster,margin=F,maxpixels=res,ylab="",xlab="",useRaster=T,scales=list(cex=1,y=list(at=c(2,6,10))),asp=1)+
  layer(sp.lines(coast,col="white",lwd=.5),under=F)


## draw it
png(width=2400,height=1600,res=200,pointsize=12,type="cairo-png",file="manuscript/figures/Seasonality.png")
print(k1,position=c(0,.5,.35,1),more=T) #legend
print(g1,position=c(0.35,.5,1,1),more=T) #global
print(r1,position=c(.46,0,1,.5),more=T)
print(r2,position=c(0,0,.48,.5),more=F)
## add panel labels
gp=gpar(fontsize=18, col="black",fontface="bold")
grid.text("a",x=0.05,y=.995,just=c("left","top"),gp=gp)
grid.text("b",x=0.35,y=.995,just=c("left","top"),gp=gp)
grid.text("c",x=0.05,y=.5,just=c("left","top"),gp=gp)
grid.text("d",x=0.51,y=.5,just=c("left","top"),gp=gp)
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
res=1e7

r="Venezuela2"
# ylab.right = "Cloud Frequency (%)",par.settings = list(layout.widths = list(axis.key.padding = 0.1,axis.left=0.6,ylab.right = 3,right.padding=2)),
pars=list(layout.heights=list(key.bottom=2,key.top=1),layout.widths = list(axis.key.padding = 3,axis.left=1))
cf_r=crop(cf_mean,regs[[r]])
p1=levelplot(cf_r,col.regions=colR(n),at=seq(0,100,len=99),
             ylab="",xlab="",
    colorkey=F,
    cuts=99,margin=F,maxpixels=1e7,par.settings = pars)+
  layer(sp.lines(coast,col="black",lwd=.5),under=F)

p2=levelplot(crop(gewex,regs[[r]]),col.regions=colR(n),at=seq(0,1,len=99),cuts=99,margin=F,max.pixels=1e6,
             ylab="",xlab="",
             colorkey=list(space="top",width=1,height=.75,labels=list(labels=c(0,50,100),at=c(.0,.5,1))),
    par.settings = pars)+
  layer(sp.lines(coast,col="black",lwd=.5),under=F)

tmap=crop(wc_map,regs[[r]])
p3=levelplot(tmap,col.regions=rev(terrain.colors(n)),cuts=100,at=seq(tmap@data@min,tmap@data@max,len=100),margin=F,maxpixels=1e6,
    colorkey=list(space="bottom",height=.75,width=1),xlab="",ylab="",main=names(regs)[r],useRaster=T,
    par.settings = pars)+
  layer(sp.lines(coast,col="black",lwd=.5),under=F)

p4=levelplot(crop(wc_dem,regs[[r]]),col.regions=terrain.colors(n),cuts=99,margin=F,max.pixels=1e6,
             ylab="",xlab="",maxpixels=1e7,
             colorkey=list(space="bottom",height=.75,width=1,labels=list(labels=c(0,2500,5000),at=c(0,2500,5000))),
    par.settings = pars,scales=list(y=list(at=c(2,6,10))))+
  layer(sp.lines(coast,col="black",lwd=.5),under=F)


png("manuscript/figures/Resolution.png",width=1000,height=2000,res=300,pointsize=47,bg="white")

trellis.par.set(my.theme)
#pdf("output/mod09_resolution.pdf",width=11,height=8.5)
#print(c("a    MODCF (%)"=p1,"b     PATMOS-x GEWEX (%)"=p2,"c    WorldClim Precip (mm)"=p3,"d    Elevation (m)"=p4,x.same=T,y.same=T,merge.legends=T,layout=c(2,2)))
print(c("a     PATMOS-x GEWEX (%)"=p2,"b    MODCF (%)"=p1,"c    Elevation (m)"=p4,x.same=T,y.same=T,merge.legends=T,layout=c(1,3)))
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
v1=xyplot(lat~lon,groups=era, data=st@data,asp=.75,auto.key=F,
          ylab="Latitude",xlab="Longitude",
          par.settings=list(superpose.symbol =list(cex=.5,pch=16,col=c("black","red"))))+
  layer(sp.lines(coast,col="black",lwd=1),under=T)

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
