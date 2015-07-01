
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
      outfile=paste0("data/MCD09_deriv/",sub(".tif","_Q20_hotspots.tif",basename(qprods[[i]])))     
      
      if(file.exists(outfile)) return(NA)
      tq=qsr[qsr$name==basename(qprods[[i]]),]#sub(".vrt",".tif",basename(tfile)),]        
      ## crop to region of interest (removing antarctica)
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
hinter=raster("data/MCD09_deriv//inter_Q20_hotspots.tif")
NAvalue(hinter)=255
hintra=raster("data/MCD09_deriv//intra_Q20_hotspots.tif")
NAvalue(hintra)=255
gain(hintra)=10
hmean=raster("data/MCD09_deriv//meanannual_Q20_hotspots.tif")
NAvalue(hmean)=255
gain(hmean)=100
  
hotspot=hinter+hintra+hmean
writeRaster(hotspot,file="data/MCD09_deriv/combined_Q20_hotspots.tif",overwrite=T,
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
  d2=d1 %.% (function(x) filter(x, complete.cases(x)))()
  
  save(d2,file="output/allvals.Rdata")
}

## Table of unique class combinations used in setting up color scheme for plots below
uvals=data.frame(table(values(hotspot)))
colnames(uvals)=c("class","count")

uvals$mean=substr(uvals$class,1,1)
uvals$intra=substr(uvals$class,2,2)
uvals$inter=substr(uvals$class,3,3)
write.csv(uvals,"output/hotspot_Q20_colors.csv",row.names=F)


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
uvals$col[uvals$name3==11]="green"  #0.2
## reset low and medium mean cloud to grey
#uvals$col[uvals$mean<3]=grey(0.8)  #0.2


## transparency to map
tcol=rgb2hsv(col2rgb(uvals$col))
tcol["s",]=tcol["s",]*uvals$mean/4
uvals$col=hsv(tcol["h",],tcol["s",],tcol["v",])


uvals$name=paste(uvals$ninter,uvals$nintra,uvals$nmean,sep=" ")
uvals$name2=paste(ifelse(uvals$ninter=="","Medium",as.character(uvals$ninter)),
                  ifelse(uvals$nintra=="","Medium",as.character(uvals$nintra)),
                  ifelse(uvals$nmean=="","Medium",as.character(uvals$nmean)),sep=" ")

uvals$include=apply(uvals[,c("ninter","nintra","nmean")],1,function(x) !any(x==""))
uvals$classid=1:nrow(uvals)

## build color table
cols=rep(grey(0.5),2^16)  #0.2
cols[uvals$class+1]=uvals$col #add 1 to account for 0-base
cols[0]=rgb(1,1,1,0)
hotspot@legend@colortable=cols


## add color table to hotspot via VRT
# create a copy
hotspot2=hotspot

## build 8-bit colrs
cols8=rep(grey(0.5),256)  #0.2
cols8[uvals$classid+1]=uvals$col #add 1 to account for 0-base
cols8[0]=rgb(1,1,1,0)
hotspot2@legend@colortable=cols8

# udpate with id number
values(hotspot2)=uvals$classid[match(values(hotspot),uvals$class)]
# write out a 8-bit version
writeRaster(hotspot2,file="data/MCD09_deriv/combined_Q20_hotspots2.tif",overwrite=T,
            datatype="INT1U",options=c("COMPRESS=LZW","PREDICTOR=2"))
# add color table using cols above
writeRasterCT(raster="data/MCD09_deriv/combined_Q20_hotspots2.tif",cols8,outfile="data/MCD09_deriv/combined_Q20_hotspots.vrt")
# convert to 3-band RGB 8-bit geotif
system(paste("gdal_translate -a_nodata 0 -co COMPRESS=LZW -co PREDICTOR=2 -of GTIFF -expand rgb  -ot Byte ",
  " data/MCD09_deriv/combined_Q20_hotspots.vrt ",
  " data/MCD09_deriv/combined_Q20_hotspots_RGB.tif")) 


#hotspot_vrt=raster("data/MCD09_deriv/combined_Q20_hotspots.vrt")
#plot(hotspot_vrt)
#NAvalue(hotspot2)=0





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
  gplot(hotspot,maxpixels=2e6)+geom_raster(aes(fill=value))+
  scale_fill_gradientn(colours=cols,guide=F,na.value="white",
                       values = 0:(length(cols)-1), 
                       rescaler = function(x, ...) x, oob = identity)+
  coord_equal()+ylim(c(-60,89))+ylab("")+xlab("")+blanktheme

#ggpairs(d1[1:3])
d3=d2[sample(1:nrow(d2),10000000),]

ginter=
  ggplot(d3,aes(y=inter/100,x=mean/100))+
  geom_point(aes(colour=uvals$col[match(class,uvals$class)],
                 order = uvals$mean[match(class,uvals$class)]),size=.3)+
  scale_y_continuous(breaks=c(0,20,40),limits=c(0,40),labels=NULL)+
  scale_x_continuous(breaks=c(0,50,100))+
  scale_colour_identity()+
  guides(colour=FALSE)+
  # stat_bin2d()+coord_equal()+
  geom_hline(yintercept=qsr[2,"Q0.1"]/100,size=.2)+
#  geom_hline(yintercept=qsr[2,"Q0.9"]/100,size=.2)+
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
  scale_x_continuous(breaks=c(0,50,100))+
  guides(colour=FALSE)+
  # stat_bin2d()+coord_equal()+
  geom_hline(yintercept=qsr[3,"Q0.1"]/100,size=.2)+
# geom_hline(yintercept=qsr[3,"Q0.9"]/100,size=.2)+
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
  scale_x_continuous(breaks=c(0,15,30))+
  guides(colour=FALSE)+
  # stat_bin2d()+coord_equal()+
  geom_hline(yintercept=qsr[3,"Q0.1"]/100,size=.2)+
#  geom_hline(yintercept=qsr[3,"Q0.9"]/100,size=.2)+
  geom_vline(xintercept=qsr[2,"Q0.1"]/100,size=.2)+
#  geom_vline(xintercept=qsr[2,"Q0.9"]/100,size=.2)+
  labs(y="Intra-annual",
       x="Inter-annual")+
  theme(strip.text = element_text(size=8),
        plot.title=element_text(size=12),
        panel.background=element_blank(),
        panel.border=element_rect(colour="black",fill="transparent"),
        plot.background=element_rect(fill="white"),
        plot.margin=unit(c(0,0,0,0), "npc"))


## create coastal plotting object
p_coast=  geom_polygon(data=gcoast,aes(x=long,y=lat,group=group),col="black",fill="transparent",size=.05)
p_coast2=  geom_polygon(data=gcoast2,aes(x=long,y=lat,group=group),col=grey(.2),fill="transparent",size=.05)

hotspot_region="EastAfrica"
#ggplot()+p_coast2+
#  coord_cartesian(xlim=c(regs[[hotspot_region]]@xmin,regs[[hotspot_region]]@xmax),
#                  ylim=c(regs[[hotspot_region]]@ymin,regs[[hotspot_region]]@ymax))
 

  ## regional box to illustrate location of inset
p_reg=geom_rect(xmax=regs[[hotspot_region]]@xmax, 
                xmin=regs[[hotspot_region]]@xmin,
                ymax=regs[[hotspot_region]]@ymax,
                ymin=regs[[hotspot_region]]@ymin,col="red",fill="transparent",size=.15)



## global plot with coast
p_hot_coast= p_hot+p_coast+p_reg


## Regional plot 
r_hotspot=crop(hotspot,regs[[hotspot_region]])

## load kenya EBAs
eba=readOGR("/media/data/Cloud/data/src/eba/ke_eba","ke_eba")

p_hot_reg=
  gplot(r_hotspot,maxpixels=2e6)+
  geom_raster(aes(fill=value))+
  scale_fill_gradientn(colours=cols,guide=F,na.value="white",
                       values = 0:(length(cols)-1), 
                       rescaler = function(x, ...) x, oob = identity)+
  ylab("")+xlab("")+blanktheme+p_coast2+
  coord_cartesian(xlim=c(regs[[hotspot_region]]@xmin,regs[[hotspot_region]]@xmax),
                  ylim=c(regs[[hotspot_region]]@ymin,regs[[hotspot_region]]@ymax))+
  theme(panel.border=element_rect(colour="red",fill="transparent"))+
  geom_path(aes(x=long,y=lat,group=group,order=order),data=fortify(eba),fill="transparent",size=.2)


## Draw it!
png("manuscript/figures/Q20_Hotspots6.png",width=3500,height=2000,
    res=600,bg="white")
print(p_hot_coast,vp=viewport(x=0,y=0.28,width=1,height=.7,just=c("left","bottom")))
print(p_key,vp=viewport(x=0,y=0,width=0.32,height=.3,just=c("left","bottom")))
print(ginter,vp=viewport(x=0.33,y=0,width=0.32,height=.3,just=c("left","bottom")))
print(gintra,vp=viewport(x=0.66,y=0,width=0.32,height=.3,just=c("left","bottom")))
## print region blowup
print(p_hot_reg,vp=viewport(x=0.08,y=.29,width=0.2,height=.35,just=c("left","bottom")))
## panel labels
pushViewport(viewport())
tgp=gpar(cex = 1, col = "black")
grid.text(label = "a" ,x = 0.08,y = .98,gp = tgp)
grid.text(label = "b" ,x = 0.09,y = .61,gp = tgp)
grid.text(label = "c" ,x = 0.09,y = .33,gp = tgp)
grid.text(label = "d" ,x = 0.39,y = .33,gp = tgp)
grid.text(label = "e" ,x = 0.73,y = .33,gp = tgp)
dev.off()


### Little color figure for earth engine
ggplot(filter(uvals,include)

       
d4=d2[sample(1:nrow(d2),1000000),]

       p_littlekey=       ggplot(d4,aes(x=inter/100,y=intra/100))+
         geom_point(aes(colour=uvals$col[match(class,uvals$class)],
                        order = uvals$mean[match(class,uvals$class)]),size=.8)+
         scale_colour_identity() +
         scale_y_continuous(breaks=c(0,20,40),limits=c(0,40))+
         scale_x_continuous(breaks=c(0,15,30))+
         guides(colour=FALSE)+
         geom_hline(yintercept=qsr[3,"Q0.1"]/100,size=.2)+
         geom_vline(xintercept=qsr[2,"Q0.1"]/100,size=.2)+
         labs(y="Intra-annual",
              x="Inter-annual")+
         theme(strip.text = element_text(size=8),
               plot.title=element_text(size=12),
               panel.background=element_blank(),
               panel.border=element_rect(colour="black",fill="transparent"),
               plot.background=element_rect(fill="white"),
               plot.margin=unit(c(0,0,0,0), "npc"))
       
png("manuscript/figures/HotspotKey.png",width=250,height=250, res=150,bg="white")
  print(p_littlekey)
       dev.off()
       

######################
### Scatterplot matrix of variables
png("manuscript/figures/scatterplot.png",width=3500,height=2500,res=600,bg="white")
rasterVis::splom(allvars)
dev.off()



