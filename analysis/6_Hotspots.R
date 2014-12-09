
# Load the libraries and set working directory
source("analysis/setup.R")


## list of products to process
qprods=list("data/MCD09_deriv/meanannual.tif",
            "data/MCD09_deriv/inter.tif",
            "data/MCD09_deriv/intra.tif",
            "data/MCD09_deriv/mean_1deg_sd.tif")


####################################################################
###  Map of hotspots
  registerDoMC(12)
  qs=c(0,0.01,0.025,0.05,0.1,0.25,0.75,0.9,0.95,.975,0.99,1)
  
  ## get elevation data for land mask
  land=raster("/mnt/data/jetzlab/Data/environ/global/worldclim/alt.bil")
  NAvalue(land)=-9999
  land=crop(land,extent(greg$xlim,greg$ylim))

  system("gdal_translate -co COMPRESS=LZW /mnt/data/jetzlab/Data/environ/global/worldclim/alt.bil data/worldclim_alt.tif")
  
  ## make table of quantiles
  if(!file.exists("output/derivedquantiles.csv")){
    qsr=foreach(i=1:length(qprods),.combine=rbind.data.frame)%dopar%{  
      tfile=paste0(tempdir(),"/",sub(".tif",".vrt",basename(qprods[[i]])))     
      ## crop to region of interest (removing antarctica)
      system(paste("gdalwarp -of vrt -te ",
                   paste(c(greg$xlim[1],greg$ylim[1],greg$xlim[2],greg$ylim[2]),collapse=" "),
                   qprods[[i]],tfile))
      ## crop to land
      trast=raster(tfile)
      trast=mask(trast,land)
      NAvalue(trast)=0
      tq=quantile(trast, probs=qs,ncells=NULL,na.rm=T,names=T)
      rmr(trast)
      file.remove(tfile)
      ## calculate quantiles
      cbind.data.frame(name=basename(qprods[[i]]),t(tq))
    }
    colnames(qsr)=c("name",paste0("Q",qs))
    
  write.csv(qsr,"output/derivedquantiles.csv",row.names=F)
  }

  ## read in quantile tables
  qsr=read.csv("output/derivedquantiles.csv")
  
  qsrc=data.frame(id=  c(0, 1,  2,   3,       4:255),
                  R=   c(0, 0,  192, 255, rep(0,252)),
                  G=   c(0, 0,  192, 0  , rep(0,252)),
                  B=   c(0, 255,192, 0  , rep(0,252)),
                  ALFA=c(0, 255,255, 255, rep(0,252)))
  write.table(qsrc,"output/derivedquantiles_colors.csv",quote=F,col.names=F,row.names=F)
  
  
  
  ## reclass each raster to low-medium-high
    foreach(i=1:length(qprods))%dopar%{  
#      tfile=paste0(tempdir(),"/",sub(".tif",".vrt",basename(qprods[[i]])))
      outfile=paste0("data/MCD09_deriv/",sub(".tif","_Q10_hotspots.tif",basename(qprods[[i]])))     
      
      if(file.exists(outfile)) return(NA)
      tq=qsr[qsr$name==basename(qprods[[i]]),]#sub(".vrt",".tif",basename(tfile)),]        
      ## crop to region of interest (removing antarctica)
#      system(paste("gdalwarp -of vrt -te ",
#                   paste(c(greg$xlim[1],greg$ylim[1],greg$xlim[2],greg$ylim[2]),collapse=" "),
#                   qprods[[i]],tfile))
      system(paste0("pksetmask -ot Byte -i ",qprods[[i]],
                    " -m data/worldclim_alt.tif --operator='<' --msknodata -1000 --nodata 0 ",                    
                    " -m ",qprods[[i]]," --operator='>' --msknodata ",tq$Q1," --nodata 0 ", # missing data
                    " -m ",qprods[[i]]," --operator='<' --msknodata ",tq$Q0.05," --nodata 1 ", # low 
                    " -m ",qprods[[i]]," --operator='<' --msknodata ",tq$Q0.95," --nodata 2 ", # middle 
                    " -m ",qprods[[i]]," --operator='>' --msknodata ",tq$Q0.95," --nodata 3 ", # top 
                    " -ct output/derivedquantiles_colors.csv ",
                    " -o ",outfile))
    }
    
  }
  
  
 
## find which combinations actually exist
hinter=raster("data/MCD09_deriv//inter_Q10_hotspots.tif")
NAvalue(hinter)=255
hintra=raster("data/MCD09_deriv//intra_Q10_hotspots.tif")
NAvalue(hintra)=255
gain(hintra)=10
hmean=raster("data/MCD09_deriv//meanannual_Q10_hotspots.tif")
NAvalue(hmean)=255
gain(hmean)=100
  
hotspot=hinter+hintra+hmean
writeRaster(hotspot,file="data/MCD09_deriv/combined_Q10_hotspots.tif",overwrite=T,
            datatype="INT2U",options=c("COMPRESS=LZW","PREDICTOR=2"))


## Make table of all values for four maps to faciliate scatterplots below.
if(!file.exists("output/allvals.Rdata")){
  ## get all values
  d1=data.frame(
    mean=values(raster(qprods[1][[1]])),
    inter=values(raster(qprods[2][[1]])),
    intra=values(raster(qprods[3][[1]])),
    class=values(hotspot))
  
  ## combine to matrix and drop NAs
  d1=values(dr)
  d2=d1 %.% (function(x) filter(x, complete.cases(x)))()
  
  save(d2,file="output/allvals.Rdata")
}

## Table of unique class combinations used in setting up color scheme for plots below
uvals=data.frame(table(values(hotspot)))
colnames(uvals)=c("class","count")

uvals$mean=substr(uvals$class,1,1)
uvals$intra=substr(uvals$class,2,2)
uvals$inter=substr(uvals$class,3,3)
write.csv(uvals,"output/hotspot_Q10_colors.csv",row.names=F)


#############################################################
### Read it back in

hotspot=raster("data/MCD09_deriv/combined_Q20_hotspots.tif")
NAvalue(hotspot)=0

load("output/allvals.Rdata")

## read in quantile tables and pixels values for scatterplots
qsr=read.csv("output/derivedquantiles.csv")



## read it back in
uvals=read.csv("output/hotspot_Q20_colors.csv")
uvals$twos=apply(uvals[,c("inter","intra","mean")],1,function(x)sum(x!=2))
uvals$ninter=factor(uvals$inter,labels=c("Low","","High"))
uvals$nintra=factor(uvals$intra,labels=c("Low","","High"))
uvals$nmean=factor(uvals$mean,labels=c("Low","Medium","High"))



## Develop colors based on inter, intra, and mean annual
# uvals[,c("h","s","v")]=t(rgb2hsv(r=uvals$mean/1.5,b=uvals$inter,g=uvals$intra,max=3))
# uvals$s=uvals$twos/3
# uvals$v=uvals$mean/3
# uvals$col=hsv(h=uvals$h,s=uvals$s,v=uvals$v)
# uvals$col[uvals$class==122]=grey(0.85)  #0.2
# uvals$col[uvals$class==222]=grey(0.9)
# uvals$col[uvals$class==322]=grey(0.95)  #0.65

## Develop colors based only in inter and intra
#uvals[,c("h","s","v")]=t(rgb2hsv(r=uvals$inter,g=2,b=uvals$intra,max=3))
#uvals$col=hsv(h=uvals$h,s=uvals$s,v=uvals$v)
uvals$name3=paste(uvals$inter,uvals$intra,sep="")
uvals$col=grey(0.95)  #0.2
#uvals$col[uvals$name3==33]="red"  #0.2
#uvals$col[uvals$name3==23]="orange"  #0.2
uvals$col[uvals$name3==13]="red"  #0.2
uvals$col[uvals$name3==12]="red"  #0.2
#uvals$col[uvals$name3==32]="purple"  #0.2
uvals$col[uvals$name3==31]="blue"  #0.2
uvals$col[uvals$name3==21]="blue"  #0.2
uvals$col[uvals$name3==11]="purple"  #0.2
## reset low and medium mean cloud to grey
#uvals$col[uvals$mean<3]=grey(0.8)  #0.2


## transparency to map
tcol=rgb2hsv(col2rgb(uvals$col))
tcol["s",]=tcol["s",]*uvals$mean/4
uvals$col=hsv(tcol["h",],tcol["s",],tcol["v",])

#uvals$col=rgb(r=uvals$inter,g=uvals$intra,b=2,max=3)

#uvals$col[uvals$class==122]=grey(0.9)  #0.2
#uvals$col[uvals$class==222]=grey(0.9)
#uvals$col[uvals$class==322]=grey(0.9)  #0.65


uvals$name=paste(uvals$ninter,uvals$nintra,uvals$nmean,sep=" ")
uvals$name2=paste(ifelse(uvals$ninter=="","Medium",as.character(uvals$ninter)),
                  ifelse(uvals$nintra=="","Medium",as.character(uvals$nintra)),
                  ifelse(uvals$nmean=="","Medium",as.character(uvals$nmean)),sep=" ")

uvals$include=apply(uvals[,c("ninter","nintra","nmean")],1,function(x) !any(x==""))

## build color table
cols=rep(grey(0.5),2^16)  #0.2
cols[uvals$class+1]=uvals$col #add 1 to account for 0-base
cols[0]=rgb(1,1,1,0)
hotspot@legend@colortable=cols

## hotspot map

tsize=20

blanktheme=theme(legend.position="none")+theme(axis.line=element_blank(),
                                               axis.text.x=element_blank(),
                                               axis.text.y=element_blank(),
                                               axis.ticks=element_blank(),
                                               axis.title.x=element_blank(),
                                               axis.title.y=element_blank(),
                                               #legend.position="none",
                                               panel.background = element_rect(fill = 'white', colour = 'white'),
                                               panel.border=element_rect(fill="transparent",linetype = "solid", colour = "transparent"),
                                               panel.grid.major=element_blank(),
                                               panel.grid.minor=element_blank(),
                                               plot.background=element_blank(),
                                               plot.margin=unit(c(0,0,0,0), "npc"),
                                               text=element_text(size=tsize),
                                               strip.background = element_blank(),
                                               strip.text = element_blank(),
                                               text = element_text(size=tsize))

p_hot=
  gplot(hotspot,maxpixels=2e4)+geom_raster(aes(fill=value))+
  scale_fill_gradientn(colours=cols,guide=F,na.value="white",
                       values = 0:(length(cols)-1), 
                       rescaler = function(x, ...) x, oob = identity)+
  coord_equal()+ylim(c(-60,89))+ylab("")+xlab("")+blanktheme

# p_key=  
#   ggplot(uvals,
#        aes(x=ninter,y=nintra,colour=col))+
# #  geom_raster()+
#   geom_point(aes(size=count))+
#   scale_size_area(trans="log",max_size=10,breaks=c(0,1e2,1e4,1e6,1e8),name="Count")+
#   scale_color_identity()+
#   facet_grid(~nmean)+scale_fill_identity() +
#   labs(title = "Mean Annual Cloud Cover", 
#        y="Intra-annual\nVariability",
#       x="Inter-annual Variability")+
#   theme(strip.text = element_text(size=8),
#         plot.title=element_text(size=12),
#         plot.background=element_blank(),
#         plot.margin=unit(c(0,.10,0,0), "npc"))

# p_key=  
#   ggplot(uvals,
#          aes(x=ninter,y=nintra,fill=col))+
#     geom_raster()+coord_equal()+
#   facet_grid(~nmean)+scale_fill_identity() +
#   labs(title = "Mean Annual Cloud Cover", 
#        y="Intra-annual\nVariability",
#        x="Inter-annual Variability")+
#   theme(strip.text = element_text(size=8),
#         plot.title=element_text(size=12),
#         plot.background=element_blank(),
#         plot.margin=unit(c(0,.10,0,0), "npc"))
# 
#  ggplot(uvals,
#         aes(x=ninter,y=nintra,colour=col))+
#    geom_point(pch=15,size=5,guide="legend")+
#    scale_colour_identity(name="Low                        Medium                       High",
#                          labels=uvals$name2,guide="legend")+
#    guides(col = guide_legend(nrow = 4))

#ggpairs(d1[1:3])
d3=d2[sample(1:nrow(d2),100000),]

ginter=
  ggplot(d3,aes(y=inter/100,x=mean/100))+
  geom_point(aes(colour=uvals$col[match(class,uvals$class)],
                 order = uvals$mean[match(class,uvals$class)]),size=.3)+
  scale_y_continuous(breaks=c(0,20,40),limits=c(0,40),labels=NULL)+
  scale_colour_identity()+
  guides(colour=FALSE)+
  # stat_bin2d()+coord_equal()+
  geom_hline(yintercept=qsr[2,"Q0.1"]/100,size=.2)+
  geom_hline(yintercept=qsr[2,"Q0.9"]/100,size=.2)+
  geom_vline(xintercept=qsr[1,"Q0.1"]/100,size=.2)+
  geom_vline(xintercept=qsr[1,"Q0.9"]/100,size=.2)+
  labs(y="Inter-annual",
       x="Mean")+
  theme(strip.text = element_text(size=8),
        plot.title=element_text(size=12),
        panel.background=element_blank(),
        panel.border=element_rect(colour="black",fill="transparent"),
        plot.background=element_rect(fill="white"),
        plot.margin=unit(c(0,0,0,0), "npc"))

gintra=
  ggplot(d3,aes(y=intra/100,x=mean/100))+
  geom_point(aes(colour=uvals$col[match(class,uvals$class)],
                 order = uvals$mean[match(class,uvals$class)]),size=.3)+
  scale_colour_identity() +
  scale_y_continuous(breaks=c(0,20,40),limits=c(0,40),labels=NULL)+
  guides(colour=FALSE)+
  # stat_bin2d()+coord_equal()+
  geom_hline(yintercept=qsr[3,"Q0.1"]/100,size=.2)+
  geom_hline(yintercept=qsr[3,"Q0.9"]/100,size=.2)+
  geom_vline(xintercept=qsr[1,"Q0.1"]/100,size=.2)+
  geom_vline(xintercept=qsr[1,"Q0.9"]/100,size=.2)+
  labs(y="Intra-annual",
       x="Mean")+
  theme(strip.text = element_text(size=8),
        plot.title=element_text(size=12),
        panel.background=element_blank(),
        panel.border=element_rect(colour="black",fill="transparent"),
        plot.background=element_rect(fill="white"),
        plot.margin=unit(c(0,0,0,0), "npc"))


#d2=data.frame(d1[sample(1:nrow(d1),100000),])

p_key=
  ggplot(d3,aes(x=inter/100,y=intra/100))+
  geom_point(aes(colour=uvals$col[match(class,uvals$class)],
                 order = uvals$mean[match(class,uvals$class)]),size=.3)+
  scale_colour_identity() +
  scale_y_continuous(breaks=c(0,20,40),limits=c(0,40))+
    guides(colour=FALSE)+
  # stat_bin2d()+coord_equal()+
  geom_hline(yintercept=qsr[3,"Q0.1"]/100,size=.2)+
  geom_hline(yintercept=qsr[3,"Q0.9"]/100,size=.2)+
  geom_vline(xintercept=qsr[2,"Q0.1"]/100,size=.2)+
  geom_vline(xintercept=qsr[2,"Q0.9"]/100,size=.2)+
  labs(y="Intra-annual",
       x="Inter-annual")+
  theme(strip.text = element_text(size=8),
        plot.title=element_text(size=12),
        panel.background=element_blank(),
        panel.border=element_rect(colour="black",fill="transparent"),
        plot.background=element_rect(fill="white"),
        plot.margin=unit(c(0,0,0,0), "npc"))



fcoast <- fortify(map_data("world"))#, plot = FALSE, fill = TRUE)

p_hot_coast=
  p_hot+
  geom_polygon(data=fcoast,aes(x=long,y=lat,group=group),col="black",fill="transparent",size=.05)

png("manuscript/figures/Q20_Hotspots3.png",width=3500,height=2000,
    res=600,bg="white")
#print(p_hot_coast,vp=viewport(x=0,y=.4,width=1,height=.6,just=c("left","bottom")))
#print(p_key,vp=viewport(x=0,y=0,width=1,height=.45,just=c("left","bottom")))
print(p_hot_coast,vp=viewport(x=0.05,y=0,width=1,height=1,just=c("left","bottom")))
#print(p_key,vp=viewport(x=0,y=0,width=0.3,height=.56,just=c("left","bottom")))
#grid.arrange(p_key, ginter,gintra, ncol=3),
print(p_key,vp=viewport(x=0,y=0,width=0.33,height=.3,just=c("left","bottom")))
print(ginter,vp=viewport(x=0.33,y=0,width=0.33,height=.3,just=c("left","bottom")))
print(gintra,vp=viewport(x=0.66,y=0,width=0.33,height=.3,just=c("left","bottom")))

## panel labels
pushViewport(viewport())
tgp=gpar(cex = 1, col = "black")
grid.text(label = "a" ,x = 0.01,y = .98,gp = tgp)
grid.text(label = "b" ,x = 0.01,y = .62,gp = tgp)
dev.off()



######################
### Scatterplot matrix of variables
png("manuscript/figures/scatterplot.png",width=3500,height=2500,res=600,bg="white")
rasterVis::splom(allvars)
dev.off()









##### Old junk below
# function to assign each grid cell to row in cgrid
getid=function(x) ifelse(any(is.na(x)),NA,which(apply(cgrid, 1, function(y) all(y == x))))
fgetid=function(x,...) calc(x,getid,...)

hs2=stack(crop(hs,regs[["Venezuela"]]))
NAvalue(hs2)=255

beginCluster(10)

clusterR(hs2,fun=fgetid,file="data/MCD09_deriv/combined_hotspots.tif",overwrite=T)

endCluster(10)
#levelplot(hs2)


tr=raster("data/MCD09_deriv/combined_hotspots.tif")
plot(tr)

table(values(tr))

}


### Build the tiles to process
jobs=tilebuilder(xmin=-180,xmax=180,ymin=-90,ymax=90,size=30,overlap=0)


file="data/MCD09_deriv/combined_hotspots.tif"
if(!file.exists(paste0(datadir,"/mcd09focal"))) dir.create(paste0(datadir,"/mcd09focal"))

registerDoMC(20)

foreach( i=1:nrow(jobs), .options.multicore=list(preschedule=FALSE)) %dopar% {
  
  toutfile1=paste(datadir,"/mcd09focal/", sub(".tif","",basename(file)),"_",jobs$tile[i],"_region.vrt",sep="")
  toutfile3=paste(datadir,"/mcd09focal/", sub(".tif","",basename(file)),"_",jobs$tile[i],"_.tif",sep="")
  
  if(file.exists(toutfile2)) {writeLines(paste(toutfile,"Exists, moving on"));return(NULL)}
  writeLines(paste("Starting: ",basename(toutfile1)," tile:",jobs$tile[i]," ( ",i," out of ",nrow(jobs),")"))
  ## crop to buffered region
  system(paste("gdalwarp -of vrt -ot Int16 -srcnodata 65535 -dstnodata -32768 -te ",
               paste(select(jobs,xminb,yminb,xmaxb,ymaxb)[i,],collapse=" "),file,toutfile1))
  calc(hs,fun=getid,file="data/MCD09_deriv/combined_hotspots.tif",overwrite=T)
  
  ## remove temporary files
  file.remove(toutfile1,toutfile2)
  print(paste("Finished Temporary File: ",basename(toutfile1)))
}


system(paste0("gdalwarp ",datadir,"/mcd09focal/*_sdcrop.tif data/MCD09_deriv/mean_1deg_sd_uncompressed.tif"))
system(paste("gdal_translate -co COMPRESS=DEFLATE -co ZLEVEL=9 -co BIGTIFF=YES -co COMPRESS=LZW -co PREDICTOR=2",
             " data/MCD09_deriv/mean_1deg_sd_uncompressed.tif data/MCD09_deriv/mean_1deg_sd.tif"))

