### Figures and tables for MOD09 Cloud Manuscript
source("analysis/setup.R")

### Load data
cf_mean=raster("data/MCD09_deriv/MCD09_meanannual.tif")

## Figures
n=100
at=seq(0,100,length=n)
colr=colorRampPalette(c("black","green","red"))
cols=colr(n)


res=1e5
greg=list(ylim=c(-60,84),xlim=c(-180,180))


### load bias correction data
## uncorrected data
mod=raster(list.files("/mnt/data2/projects/cloud/mcd09tif",pattern="MOD09_01_mean.tif",full=T))
modc=raster(list.files("/mnt/data2/projects/cloud/mcd09ctif",pattern="MOD09_01_mean.tif",full=T))
modc3=raster(list.files("/mnt/data2/projects/cloud/mcd09ctif",pattern="MOD09_03_mean.tif",full=T))

## interannual variability
inter=raster("data/MCD09_deriv/inter_old.tif")
gain(inter)=0.01
intra=raster("data/MCD09_deriv/intra.tif")
gain(intra)=0.01


regs=list(
  Sahara=extent(c(-5,5,5,15))
)

r="Sahara"
dt=crop(stack(mod,modc),regs[[r]])/100
names(dt)=c("A)Uncorrected_Terra ","B)Corrected_Terra")

## figure settings
pres.theme = trellis.par.get()
pres.theme$strip.background=list(col="transparent")
pres.theme$par.ylab.text$cex=2
pres.theme$par.xlab.text$cex=2
pres.theme$axis.text$cex=2
pres.theme$layout.widths$key.right=1.5

##################
### Presentation/Poster
res=1e7

pdf("poster/figures/MeanAnnualCF.pdf",width=11,height=8.5,pointsize=36)

trellis.par.set(pres.theme)
pars=list(layout.heights=list(strip=1.5,key.bottom=2,key.top=1),layout.widths = list(axis.key.padding = 3,axis.left=1))
levelplot(cf_mean,col.regions=colR(n),cuts=99,at=seq(0,100,len=100),colorkey=list(space="right",adj=1),
          main="Mean Annual Cloud Frequency (%)",
          panel=panel.levelplot.raster,margin=F,maxpixels=res,ylab="",xlab="",useRaster=T,ylim=greg$ylim)+
  layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col="white"),under=T)+
  layer(sp.lines(hcoast,col=grey(0.3),lwd=.4),under=F)
dev.off()


res=1e7

pdf("poster/figures/InterAnnual.pdf",width=11,height=8.5,pointsize=36)
trellis.par.set(pres.theme)
pars=list(layout.heights=list(strip=1.5,key.bottom=2,key.top=1),layout.widths = list(axis.key.padding = 3,axis.left=1))
levelplot(inter,col.regions=bgyrp(n),cuts=n-1,colorkey=list(space="right",adj=1),
          main="Interannual Standard Deviation",
          panel=panel.levelplot.raster,margin=F,maxpixels=res,ylab="",xlab="",useRaster=T,ylim=c(-60,69))+
  layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col="white"),under=T)+
  layer(sp.lines(hcoast,col="black",lwd=.4),under=F)
dev.off()


#system(paste("gdalinfo -stats ",intra@file@name))
pdf("poster/figures/IntraAnnual.pdf",width=11,height=8.5,pointsize=36)
trellis.par.set(pres.theme)
pars=list(layout.heights=list(strip=1.5,key.bottom=3,key.top=1),layout.widths = list(axis.key.padding = 4,axis.bottom=4,axis.left=1))
#at_intra=seq(0,155,len=n-1)
levelplot(intra,col.regions=bgyrp(n),cuts=n-1,colorkey=list(space="right",adj=1),
          main="Intraannual Standard Deviation",
          panel=panel.levelplot.raster,
          subscripts=T,margin=F,maxpixels=res,xlab="",ylab="",useRaster=T,ylim=c(-60,69))+
  layer(panel.polygon(x=c(-180,-180,180,180),y=c(-90,90,90,-90),col="white"),under=T)+
  layer(sp.lines(hcoast,col="black",lwd=.4),under=F)
dev.off()


### Example monthly mean
pdf("poster/figures/MarchMean.pdf",width=11,height=8.5,pointsize=36)
trellis.par.set(pres.theme)

pars=list(layout.heights=list(strip=1.5,key.bottom=2,key.top=1),layout.widths = list(axis.key.padding = 3,axis.left=1))
bcols=colorRampPalette(c("#08306b","#0d57a1","#2878b8","#4997c9","#72b2d7","#a2cbe2","#c7dcef","#deebf7","#f7fbff"))
p_mar=levelplot(modc3,col.regions=bcols(100),at=seq(0,100,len=101),ylim=greg$ylim,
                    panel=panel.levelplot.raster,ylab.right="MODIS Cloud Frequency (%)",
                    cuts=99,margin=F,maxpixels=1e6,par.settings = pars)+layer(sp.lines(hcoast))
p_mar$par.strip.text=list(cex = 2)
print(p_mar)
dev.off()

###############################################

pdf("poster/figures/PosterFigures.pdf",useDingbats=F,width=11,height=8.5,pointsize=36)
trellis.par.set(pres.theme)

xyplot(lat~lon,data=cldm[cldm$month==1,],pch=16,cex=.5,col="black",
              ylab="Latitude",xlab="Longitude",ylim=greg$ylim)+
layer(sp.polygons(land,fill="#c7dcef",lwd=.1),col=grey(0.5),under=T)



### heatmap of mod09 vs. NDP for all months
hmcols=colorRampPalette(c("grey","blue","red","purple"))
#hmcols=colorRampPalette(c(grey(.8),grey(.3),grey(.2)))
tr=c(0,50)
colkey <- draw.colorkey(list(col = hmcols(tr[2]), at = tr[1]:tr[2],labels=list(at=c(0,25,50),labels=c(0,25,50)),height=.4))

xyplot(cld~mod09|seas,data=cldm,panel=function(x,y,subscripts){
  n=50
  bins=seq(0,100,len=n)
  tb=melt(as.matrix(table(
    x=cut(x,bins,labels=bins[-1]),
    y=cut(y,bins,labels=bins[-1]))))
  qat=unique(tb$value)
  print(max(qat))
  qat=tr[1]:tr[2]#unique(tb$value)
  panel.levelplot.raster(tb$x,tb$y,tb$value,at=qat,col.regions=c("transparent",hmcols(length(qat))),subscripts=1:nrow(tb))
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

##########
## Bias Correction
pars=list(layout.heights=list(strip=1.5,key.bottom=3,key.top=2),layout.widths = list(axis.key.padding = 3,axis.left=1))
bcols=colorRampPalette(c("#08306b","#0d57a1","#2878b8","#4997c9","#72b2d7","#a2cbe2","#c7dcef","#deebf7","#f7fbff"))
p_bias1=levelplot(dt,ylab="", panel=panel.levelplot.raster,xlab="MODIS Cloud Frequency (%)",cuts=99,
          margin=F,col.regions=bcols(100),at=seq(0,100,len=101),colorkey=list(space="bottom"),maxpixels=5e5,par.settings = pars)#+layer(sp.lines(hcoast))
p_bias1$par.strip.text=list(cex = 2)
print(p_bias1)

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

