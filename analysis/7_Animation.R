## script to visualize cloud frequency data
source("analysis/setup.R")



#####  Daily animation
clda=stack(list.files("/mnt/data2/scratch/tmp/MCD09/"))
cols=colorRampPalette(c("#000000","#00FF00","#FF0000"))#"black","blue","red"))
foreach(i=1:12)%dopar% {
    print(i)
    r=mc[[i]]
    d1=(extent(coast2)@ymax-extent(coast2)@ymin)/(extent(coast2)@xmax-extent(coast2)@xmin)  # get y:x ratio to set output
    png(paste("output/CF_mean_",sprintf("%02d", i),".png",sep=""),width=1920,height=round(1920*d1),res=300,pointsize=46,bg="white")
    print(levelplot(r,col.regions=cols(100),at=seq(1,100,len=100),margin=F,maxpixels=2.1e6,ylim=c(extent(coast2)@ymin,extent(coast2)@ymax),
                    main=paste(month.name[i]),cex.main=3,scales=list(draw=F),cuts=99,ylab="",xlab="")+
                        layer(panel.polygon(c(extent(world)@xmin,extent(world)@xmin,extent(world)@xmax,extent(world)@xmax),y=c(extent(world)@ymin,extent(world)@ymax,extent(world)@ymax,extent(world)@ymin),col="black"),under=T)+
                        layer(sp.lines(coast2,col="black"),under=F)+
                        layer(sp.lines(world,col="black"),under=F))

    dev.off()
}




####################################################################
### Regional Comparisons
## Compare with worldclim and NPP
wc=stack(as.list(paste("/mnt/data/jetzlab/Data/environ/global/worldclim/prec_",1:12,".bil",sep="")))
#wc_map=stack(as.list(paste("/mnt/data/jetzlab/Data/environ/global/worldclim/bio_12.bil",sep="")))

regs=list(
  Cascades=extent(c(-122.8,-118,44.9,47)),
  Hawaii=extent(c(-156.5,-154,18.75,20.5)),
  Boliva=extent(c(-71,-63,-20,-15)),
  Venezuela=extent(c(-69,-59,0,7)),
  CFR=extent(c(17.75,22.5,-34.8,-32.6)),
  Madagascar=extent(c(46,52,-17,-12))
  #reg2=extent(c(-81,-70,-4,10))
  )
# convert to sinusoidal
regs2=lapply(regs,function(r){
       r2 = as(r, "SpatialPolygons")
       proj4string(r2) <- CRS(proj4string(land))
       r3=spTransform(r2,CRS(projection(mc)))
       return(r3)
   })


## read in GEWEX 1-degree data
gewex=brick("data/gewex/CA_PATMOSX_NOAA.nc",varname="a_CA")
names(gewex)=month.name#"1-degree Cloud Frequency (PATMOS-x GEWEX AVHRR)"






r="Venezuela"

## print global map with box for region
png(paste("output/CF_mean_regbox.png",sep=""),width=1920,height=round(1920*d1),res=300,pointsize=46,bg="white")
print(levelplot(mc[[1]],col.regions=cols(100),at=seq(1,100,len=100),margin=F,maxpixels=2.1e6,ylim=c(extent(coast2)@ymin,extent(coast2)@ymax),
                main=paste(month.name[i]),cex.main=3,scales=list(draw=F),cuts=99,ylab="",xlab="")+
      layer(panel.polygon(c(extent(world)@xmin,extent(world)@xmin,extent(world)@xmax,extent(world)@xmax),y=c(extent(world)@ymin,extent(world)@ymax,extent(world)@ymax,extent(world)@ymin),col="black"),under=T)+
      layer(sp.lines(coast2,col="black"),under=F)+
      layer(sp.lines(world,col="black"),under=F)+
      layer(sp.lines(as(regs2[[r]],"SpatialLines"),col="blue",lwd=2),under=F))
dev.off()

                                        # ylab.right = "Cloud Frequency (%)",par.settings = list(layout.widths = list(axis.key.padding = 0.1,axis.left=0.6,ylab.right = 3,right.padding=2)),

    ## crop the data
    tgewex=crop(gewex,regs[[r]])
    tmap=crop(wc,regs[[r]])
    tmc=projectRaster(crop(mc,regs2[[r]]),tmap)

## set plotting parameters

    range(c(cellStats(tmc,range)),cellStats(tgewex,range)*100)
range(cellStats(tmap,range))

pcols=colorRampPalette(c("#CCFFFF","#000066","#FF3300"))#"black","blue","red"))


foreach(i=1:12)%dopar% {
    print(i)
    png(paste("output/CF_mean_",r,"_",sprintf("%02d", i),".png",sep=""),width=1920,height=1080,res=300,pointsize=46,bg="white")
    pars=list(layout.heights=list(key.bottom=2,key.top=1),layout.widths = list(axis.key.padding = 4,axis.left=1))
    my.theme = trellis.par.get()
    my.theme$strip.background=list(col="transparent")
    trellis.par.set(my.theme)
    p1=levelplot(tmc[[i]],col.regions=cols(100),at=seq(20,100,len=99),
        colorkey=list(space="bottom",width=1,height=1.2,labels=list(labels=c(25,50,75,100),at=c(25,50,75,100)),just="left"),
        cuts=99,margin=F,max.pixels=1e6,par.settings = pars,main=paste(month.name[i]),cex.main=3,xlab="",ylab="",scales=list(x=list(at=c(-68,-64,-60))))
    p2=levelplot(tgewex[[i]],col.regions=cols(100),at=seq(.20,1,len=99),cuts=99,margin=F,max.pixels=1e6,
        colorkey=F,#list(space="bottom",width=.01,height=0.00000001,labels=list(labels="",at="",cex=0.0001)),
        par.settings = pars)
    p3=levelplot(tmap[[i]],col.regions=pcols(100),cuts=100,at=seq(0,650,len=100),margin=F,maxpixels=1e6,
        colorkey=list(space="bottom",height=.5,width=1,labels=list(labels=c(0,300,600),at=c(0,300,600))),xlab="",ylab="",main=names(regs)[r],useRaster=T,
        par.settings = pars)
    p_all=c("MODIS Cloud (%)"=p1,"PATMOS-x Cloud (%)"=p2,"WorldClim Precip (mm)"=p3,x.same=T,y.same=T,merge.legends=T,layout=c(3,1))
    print(p_all)
    dev.off()
}
