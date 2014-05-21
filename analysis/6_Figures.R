### Figures and tables for MOD09 Cloud Manuscript
source("analysis/setup.R")


### Load data
cf_mean=raster("data/MCD09_deriv/MCD09_meanannual.tif")

## Figures
n=100
at=seq(0,100,length=n)
colr=colorRampPalette(c("black","green","red"))
cols=colr(n)

## set plotting parameters
my.theme = trellis.par.get()
my.theme$strip.background=list(col="transparent")
trellis.par.set(my.theme)

pdf("manuscript/figures/Figures.pdf",width=11,height=8.5,pointsize=14)
#png("manuscript/figures/CF_Figures_%03d.png",width=5000,height=4000,res=600,pointsize=42,bg="white")

## set plotting parameters
my.theme = trellis.par.get()
my.theme$strip.background=list(col="transparent")
trellis.par.set(my.theme)

res=1e5
greg=list(ylim=c(-60,84),xlim=c(-180,180))
    
## Figure 1: 4-panel summaries
#- Annual average
levelplot(cf_mean,col.regions=colR(n),cuts=99,at=seq(0,100,len=100),colorkey=list(space="right",adj=1),
          panel=panel.levelplot.raster,margin=F,maxpixels=res,ylab="",xlab="",useRaster=T,ylim=greg$ylim)+
    layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col="white"),under=T)+
    layer(sp.lines(hcoast,lwd=.2,),under=F)

## Mean annual with validation stations
levelplot(cf_mean,col.regions=colr(n),cuts=99,at=seq(0,100,len=100),colorkey=list(title="Cloud Frequency (%)",space="bottom",adj=1),
  margin=F,maxpixels=res,ylab="",xlab="",useRaster=T,ylim=greg$ylim)+
    layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col="black"),under=T)+
    layer(panel.xyplot(lon,lat,pch=16,cex=.3,col="black"),data=data.frame(coordinates(st)))+
  layer(sp.lines(coast,col="black"),under=F)

## Seasonal Means
levelplot(mod09s,col.regions=colr(n),cuts=100,at=seq(0,100,len=100),colorkey=list(title="Cloud Frequency (%)", space="bottom",adj=2),
  margin=F,maxpixels=res,ylab="",xlab="",useRaster=T,ylim=greg$ylim)+
    layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col="black"),under=T)+
    layer(sp.lines(coast,col="black"),under=F)


## Monthly Means
levelplot(mod09c,col.regions=colr(n),cuts=100,at=seq(0,100,len=100),colorkey=list(title="Cloud Frequency (%)", space="bottom",adj=1), 
  margin=F,maxpixels=res,ylab="Latitude",xlab="Longitude",useRaster=T,ylim=greg$ylim)+
    layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col="black"),under=T)+
    layer(sp.lines(coast,col="black"),under=F)

#- Monthly minimum
#- Monthly maximum
#- STDEV or Min-Max
p_mac=levelplot(mod09a,col.regions=colr(n),cuts=99,at=seq(0,100,length=100),margin=F,
    ylim=greg$ylim,maxpixels=res/10,colorkey=list(title="Cloud Frequency (%)", space="top",height=.75),xlab="",ylab="",useRaster=T)+
    layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col="black"),under=T)+
    layer(sp.lines(coast,col="black"),under=F)
p_max=levelplot(mod09max,col.regions=colr(n),cuts=99,at=seq(0,100,length=100),margin=F,maxpixels=res/10,
    ylim=greg$ylim,colorkey=list(space="bottom",height=.75),useRaster=T)+
    layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col="black"),under=T)+
    layer(sp.lines(coast,col="black"),under=F)
p_intra=levelplot(mod09intra,col.regions=colr(n),cuts=99,at=seq(0,40,length=100),margin=F,maxpixels=res/100,
    ylim=greg$ylim,colorkey=list(space="bottom",height=.75),useRaster=T)+
    layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col="black"),under=T)+
    layer(sp.lines(coast,col="black"),under=F)
p_inter=levelplot(mod09inter,col.regions=colr(n),cuts=99,at=seq(0,15,length=100),margin=F,maxpixels=res/100,
    ylim=greg$ylim,colorkey=list(title="Cloud Frequency (%)", space="top",height=.75),useRaster=T)+
    layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col="black"),under=T)+
    layer(sp.lines(coast,col="black"),under=F)

p3=c("Mean Cloud Frequency (%)"=p_mac,"Max Cloud Frequency (%)"=p_max,"Interannual Variability (sd)"=p_inter,"Intraannual Variability (sd)"=p_intra,x.same=T,y.same=F,merge.legends=T,layout=c(2,2))
print(p3)

bgr=function(x,n=100,br=0,c1=c("darkblue","blue","grey"),c2=c("grey","red","purple")){
    at=unique(c(seq(min(x,na.rm=T),max(x,na.rm=T),len=n)))
    bg=colorRampPalette(c1)
    gr=colorRampPalette(c2)
    return(list(at=at,col=c(bg(sum(at<br)),gr(sum(at>=br)))))
}

###################
### Presentation/Poster
pdf("poster/figures/PosterFigures.pdf",width=11,height=8.5,pointsize=36)
res=1e7

pres.theme = trellis.par.get()
pres.theme$strip.background=list(col="transparent")
pres.theme$par.ylab.text$cex=2
pres.theme$par.xlab.text$cex=2
pres.theme$axis.text$cex=2
pres.theme$layout.widths$key.right=1.5
trellis.par.set(pres.theme)

levelplot(cf_mean,col.regions=colR(n),cuts=99,at=seq(0,100,len=100),colorkey=list(space="right",adj=1),
          panel=panel.levelplot.raster,margin=F,maxpixels=res,ylab="",xlab="",useRaster=T,ylim=greg$ylim,
          scales=)+
  layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col="white"),under=T)+
  layer(sp.lines(hcoast,lwd=.2,),under=F)

xyplot(lat~lon,data=cldm[cldm$month==1,],pch=16,cex=.25,col="red",
              ylab="Latitude",xlab="Longitude",ylim=greg$ylim)+
layer(sp.polygons(hland,fill=grey(.8),lwd=.1),col=grey(0.5),under=T)


### heatmap of mod09 vs. NDP for all months
hmcols=colorRampPalette(c("grey","blue","red","purple"))
#hmcols=colorRampPalette(c(grey(.8),grey(.3),grey(.2)))
tr=c(0,137)
colkey <- draw.colorkey(list(col = hmcols(tr[2]), at = tr[1]:tr[2],labels=list(at=c(0,60,120),labels=c(0,60,120)),height=.4))

xyplot(cld_all~mod09|seas,data=cldm,panel=function(x,y,subscripts){
  n=50
  bins=seq(0,100,len=n)
  tb=melt(as.matrix(table(
    x=cut(x,bins,labels=bins[-1]),
    y=cut(y,bins,labels=bins[-1]))))
  qat=unique(tb$value)
  print(max(qat))
  qat=tr[1]:tr[2]#unique(tb$value)
  panel.levelplot(tb$x,tb$y,tb$value,at=qat,col.regions=c("transparent",hmcols(length(qat))),subscripts=1:nrow(tb))
  #  panel.abline(0,1,col="black",lwd=2)
  panel.abline(lm(y ~ x),col="black",lwd=2)
  #  panel.ablineq(lm(y ~ x), r.sq = TRUE,at = 0.6,pos=1, offset=0,digits=2,col="blue")
  panel.text(70,10,bquote(paste(R^2,"=",.(round(summary(lm(y ~ x))$r.squared,2)))),cex=1.5)
},asp=1,scales=list(at=seq(0,100,len=3)),useRaster=T,
                    ylab.right="# of Stations",
                    ylab="Validation Mean Cloud\nAmount (%)",xlab="MODIS Cloud Frequency (%)",
                    alternating=2, strip=strip.custom(par.strip.text = list(cex = 2)),
       par.settings = list(layout.heights=list(strip=1.5)),as.table=T,
       legend= list(right = list(fun = colkey)))#+ layer(panel.abline(0,1,col="black",lwd=2))

dev.off()



bwplot(lulcc~difm,data=cldm,horiz=T,xlab="Difference (MOD09-Observed)",varwidth=T,notch=T)+layer(panel.abline(v=0))

dev.off()

####################################################################
### Regional Comparisons
## Compare with worldclim and NPP
#wc=stack(as.list(paste("/mnt/data/jetzlab/Data/environ/global/worldclim/prec_",1:12,".bil",sep="")))
wc_map=stack(as.list(paste("/mnt/data/jetzlab/Data/environ/global/worldclim/bio_12.bil",sep="")))
wc_dem=stack(as.list(paste("/mnt/data/jetzlab/Data/environ/global/worldclim/alt.bil",sep="")))


## read in GEWEX 1-degree data
gewex=mean(brick("data/gewex/CA_PATMOSX_NOAA.nc",varname="a_CA"))
names(gewex)="PATMOS-x GEWEX AVHRR"

## calculate 1-degree means of MODCF data
#MOD_gewex=gewex
#MOD_gewex@data@values=1:length(MOD_gewex@data@values)
#MOD_gewex2=zonal(mod09a,MOD_gewex,fun='mean')
png("output/Resolution_Figures_%03d.png",width=5500,height=4000,res=600,pointsize=36,bg="white")
trellis.par.set(my.theme)
#pdf("output/mod09_resolution.pdf",width=11,height=8.5)

r="Venezuela"
# ylab.right = "Cloud Frequency (%)",par.settings = list(layout.widths = list(axis.key.padding = 0.1,axis.left=0.6,ylab.right = 3,right.padding=2)),
pars=list(layout.heights=list(key.bottom=2,key.top=1),layout.widths = list(axis.key.padding = 3,axis.left=0.6))
p1=levelplot(crop(mod09a,regs[[r]]),col.regions=grey(seq(0,1,len=100)),at=seq(45,100,len=99),
    colorkey=list(space="top",width=1,height=.75,labels=list(labels=c(50,75,100),at=c(50,75,100))),
    cuts=99,margin=F,max.pixels=1e6,par.settings = pars)
p2=levelplot(crop(gewex,regs[[r]]),col.regions=grey(seq(0,1,len=100)),at=seq(.45,1,len=99),cuts=99,margin=F,max.pixels=1e6,
    colorkey=list(space="top",width=1,height=.75,labels=list(labels=c(50,75,100),at=c(.50,.75,1))),
    par.settings = pars)
tmap=crop(wc_map,regs[[r]])
p3=levelplot(tmap,col.regions=grey(seq(0,1,len=100)),cuts=100,at=seq(tmap@data@min,tmap@data@max,len=100),margin=F,maxpixels=1e6,
    colorkey=list(space="bottom",height=.75,width=1),xlab="",ylab="",main=names(regs)[r],useRaster=T,
    par.settings = pars)
p4=levelplot(crop(wc_dem,regs[[r]]),col.regions=grey(seq(0,1,len=100)),cuts=99,margin=F,max.pixels=1e6,
    colorkey=list(space="bottom",height=.75,width=1),
    par.settings = pars)#,labels=list(labels=c(1000,4000),at=c(1000,4000))))
print(c("MODCF (%)"=p1,"PATMOS-x GEWEX (%)"=p2,"WorldClim Precip (mm)"=p3,"Elevation (m)"=p4,x.same=T,y.same=T,merge.legends=T,layout=c(2,2)))


dev.off()

###  Spatial Autocorrelation





seasl=as.data.frame(rasterToPoints(seas2))
seasl$col=col$val[match(seasl$seas_vis,col$id)]

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
my.theme = trellis.par.get()
my.theme$strip.background=list(col="transparent")
trellis.par.set(my.theme)

print(l1,position=c(0,0,1,.25),more=T)
print(l2,position=c(0,.2,1,1),more=F)

dev.off()

# mangle=seq(30,360,30)*(pi/180)
plot(col$x,col$y,xaxt="n",yaxt="n",xlab="",ylab="",col=col$val,pch=16,cex=1.5,new=F,asp=1,bty="n",xlim=c(-16,16),ylim=c(-40,40))
plot(col$x[col$exists],col$y[col$exists],xaxt="n",yaxt="n",xlab="",ylab="",col=col$val[col$exists],pch=16,cex=2,new=F,asp=1,bty="n",xlim=c(-16,16),ylim=c(-45,45))
n=10; circ=seq(0,30,n) #define the circle diameters
symbols(rep(0,length(circ)),rep(0,length(circ)),circles=circ,add=T,fg="grey",lwd=2,inches=F) #draw circles
segments(0,0,30*cos(mangle),30*sin(mangle),col="grey") #draw angles
text(x=-2,y=-circ[2:n],circ[2:n],pos=4,cex=1.5) #add scale text
mon=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec") # month names
ladj=35
for(i in c(1:3,9:12)) text(x=ladj*cos(-mangle[i]),y=ladj*sin(-mangle[i]),mon[i],pos=4,cex=1.8,offset=0,srt=-mangle[i]*180/pi) #add left months
for(i in 4:8) text(x=ladj*cos(-mangle[i]),y=ladj*sin(-mangle[i]),mon[i],pos=2,cex=1.8,offset=0,srt=(-mangle[i]*180/pi)-180) # add right months


# pushViewport(viewport(x=.75,y=.65,width=.42,height=.35,name="d1"))
# grid.rect(gp=gpar(fill="white", lty=1))
# par(fig=gridFIG(),mar=c(0,0,0,0),new=T)
# popViewport()

## vector plot of seasonality
vectorplot(seas,narrows=2e3, lwd.arrows=0.6, length=unit(5e-2, 'npc'),
           maxpixels=1e5, region=TRUE, margin=FALSE,
           isField=TRUE, reverse=FALSE,
           unit='degrees', scaleSlope=TRUE,
           aspX=0.08)
