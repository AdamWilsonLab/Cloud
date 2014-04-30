### Figures and tables for MOD09 Cloud Manuscript
source("analysis/setup.R")

## calculated differences
cldm$difm=cldm$mod09-cldm$cld_all
cldm$difs=cldm$mod09sd+cldm$cldsd_all

#clda$dif=clda$mod09-clda$cld

## read in global coasts for nice plotting
library(maptools)
library(rgdal)


## Figures
n=100
at=seq(0,100,length=n)
colr=colorRampPalette(c("black","green","red"))
cols=colr(n)

## set plotting parameters
my.theme = trellis.par.get()
my.theme$strip.background=list(col="transparent")
trellis.par.set(my.theme)

#pdf("output/Figures.pdf",width=11,height=8.5)
png("manuscript/figures/CF_Figures_%03d.png",width=5000,height=4000,res=600,pointsize=42,bg="white")

## set plotting parameters
my.theme = trellis.par.get()
my.theme$strip.background=list(col="transparent")
trellis.par.set(my.theme)

res=1e6
greg=list(ylim=c(-60,84),xlim=c(-180,180))
    
## Figure 1: 4-panel summaries
#- Annual average
levelplot(mod09a,col.regions=colr(n),cuts=100,at=seq(0,100,len=100),colorkey=list(space="bottom",adj=1),
  margin=F,maxpixels=res,ylab="",xlab="",useRaster=T,ylim=greg$ylim)+
    layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col="black"),under=T)+
    layer(sp.lines(coast,col="black"),under=F)
## Mean annual with validation stations
levelplot(mod09a,col.regions=colr(n),cuts=100,at=seq(0,100,len=100),colorkey=list(title="Cloud Frequency (%)",space="bottom",adj=1),
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

cldm$resid=NA
# get residuals of simple linear model
cldm$resid[!is.na(cldm$cld_all)&!is.na(cldm$mod09)]=residuals(lm(mod09~cld_all,data=cldm))
colat=bgr(cldm$resid)
phist=histogram(cldm$resid,breaks=colat$at,border=NA,col=colat$col,xlim=c(-30,30),type="count",xlab="MODCF Residuals")#,seq(0,1,len=6),na.rm=T)
pmap=xyplot(lat~lon|month2,data=cldm,groups=cut(cldm$resid,rev(colat$at)),
       par.settings=list(superpose.symbol=list(col=colat$col)),pch=16,cex=.25,
       auto.key=F,#list(space="right",title="Difference\n(MOD09-NDP026D)",cex.title=1),asp=1,
       ylab="Latitude",xlab="Longitude")+
  layer(sp.lines(coast,col="black",lwd=.1),under=F)
print(phist,position=c(0,.75,1,1),more=T)
print(pmap,position=c(0,0,1,.78))

### heatmap of mod09 vs. NDP for all months
hmcols=colorRampPalette(c("grey","blue","red","purple"))
#hmcols=colorRampPalette(c(grey(.8),grey(.3),grey(.2)))
tr=c(0,66)
colkey <- draw.colorkey(list(col = hmcols(tr[2]), at = tr[1]:tr[2],height=.25))

xyplot(cld_all~mod09|month2,data=cldm,panel=function(x,y,subscripts){
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
  panel.text(70,10,bquote(paste(R^2,"=",.(round(summary(lm(y ~ x))$r.squared,2)))),cex=1)
},asp=1,scales=list(at=seq(0,100,len=6),useRaster=T,colorkey=list(width=.5,title="Number of Stations")),
          ylab="NDP Mean Cloud Amount (%)",xlab="MODCF Cloud Frequency (%)",
              legend= list(right = list(fun = colkey)))#+ layer(panel.abline(0,1,col="black",lwd=2))


bwplot(lulcc~difm,data=cldm,horiz=T,xlab="Difference (MOD09-Observed)",varwidth=T,notch=T)+layer(panel.abline(v=0))

dev.off()

####################################################################
### Regional Comparisons
## Compare with worldclim and NPP
#wc=stack(as.list(paste("/mnt/data/jetzlab/Data/environ/global/worldclim/prec_",1:12,".bil",sep="")))
wc_map=stack(as.list(paste("/mnt/data/jetzlab/Data/environ/global/worldclim/bio_12.bil",sep="")))
wc_dem=stack(as.list(paste("/mnt/data/jetzlab/Data/environ/global/worldclim/alt.bil",sep="")))

regs=list(
  Cascades=extent(c(-122.8,-118,44.9,47)),
  Hawaii=extent(c(-156.5,-154,18.75,20.5)),
  Boliva=extent(c(-71,-63,-20,-15)),
  Venezuela=extent(c(-69,-59,0,7)),
  CFR=extent(c(17.75,22.5,-34.8,-32.6)),
  Madagascar=extent(c(46,52,-17,-12))
  #reg2=extent(c(-81,-70,-4,10))
  )


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

