## Figures used for NASA biodiversity meeting poster

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
