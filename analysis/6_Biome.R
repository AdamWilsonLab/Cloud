### Produce summary of cloud frequency by biome

setwd("~/acrobates/adamw/projects/cloud/")


## libraries
library(rasterVis)
library(latticeExtra)
library(xtable)
library(texreg)
library(reshape)
library(caTools)
library(rgeos)
library(raster)


## Create a simplified version of the TEOW biome dataset
if(!file.exists("data/teow/biomes.shp")){
  library(geosphere)
  library(rgeos)
  teow=readOGR("/mnt/data/jetzlab/Data/environ/global/teow/official/","wwf_terr_ecos")
  biome=unionSpatialPolygons(teow,paste(teow$REALM,teow$BIOME,sep="_"), threshold=5)
  biomeid=read.csv("/mnt/data/jetzlab/Data/environ/global/teow/biome.csv",stringsAsFactors=F)
  realmid=read.csv("/mnt/data/jetzlab/Data/environ/global/teow/realm.csv",stringsAsFactors=F,na.strings = "TTTT")
  dt=data.frame(code=row.names(biome),stringsAsFactors=F)
  dt[,c("realmid","biomeid")]=do.call(rbind,strsplit(sub(" ","",dt$code),"_"))
  dt$realm=realmid$realm[match(dt$realmid,realmid$realmid)]
  dt$biome=biomeid$BiomeID[match(dt$biomeid,biomeid$Biome)]
  row.names(dt)=row.names(biome)
  biome=SpatialPolygonsDataFrame(biome,data=dt)
  ## add area and centroid to each polygon
  biome$areakm=do.call(c,mclapply(1:length(biome),function(i) {print(i); return(areaPolygon(biome[i,])/1000000)}))
  biome@data[,c("lon","lat")]=coordinates(gCentroid(biome,byid=T))
  ## add numeric biome code for rasterization
  biome=biome[order(biome$realm,as.numeric(biome$biomeid)),]
  biome$icode=1:nrow(biome)
  ## write it to disk as shapefile
  writeOGR(biome,"data/teow","biomes",driver="ESRI Shapefile",overwrite=T)
}

## rasterize the biome dataset
## create a copy with a numeric biome code
biome=readOGR("data/teow/","biomes")
bcode=unique(data.frame(icode=biome$icode,code=biome$code,realm=biome$realm,biome=biome$biome))

system(paste("gdal_rasterize -a icode -init 0 -l biomes -ot Byte -te -180 -90 180 90 -tr 0.008333333333333 -0.008333333333333",
              " -co COMPRESS=LZW -co ZLEVEL=9 -co PREDICTOR=2 ",
              " data/teow/biomes.shp  data/teow/teow.tif"))

### Summarize the cloud data by biome
foreach(m=1:12)%dopar%{
  tm=sprintf("%02d",m)
  tcloud=paste("data/MCD09/MCD09_mean_",tm,".tif",sep="")
  tbiome=paste("data/teow/teow_MCD09_mean_",tm,".tif",sep="")
  tcloudbiome=paste("data/teow/teow_MCD09_mean_",tm,".txt",sep="")
  
  ## mask biome raster using missing data in cloud dataset
  system(paste("pksetmask -i data/teow/teow.tif -m ",tcloud," -ot UInt8 ",
               "--operator='>' --msknodata 100 --nodata 0  -co COMPRESS=LZW -co PREDICTOR=2 -o ",tbiome))
  ## calculate biome-level summary metrics
  system(paste("oft-stat -i ",tcloud," -o ",tcloudbiome," -um ",tbiome," -mm"))
  ## clean up
  file.remove(tbiome)
}    

bs=do.call(rbind.data.frame,lapply(1:12,function(m){
  tm=sprintf("%02d",m)
  tcloudbiome=paste("data/teow/teow_MCD09_mean_",tm,".txt",sep="")
  td=read.table(tcloudbiome,col.names=c("icode","n","min","max","mean","sd"))  
  td$meanpsd=td$mean+td$sd
  td$meanmsd=td$mean-td$sd
  td$month=tm
  td=merge(td,bcode,by="icode")
  file.remove(tcloudbiome)
  return(td)
  }))
write.csv(bs,file="data/sum/biomesummary.csv",row.names=F)

###################################################################
### summary by biome
bs=read.csv(file="data/sum/biomesummary.csv")
bs$monthname=factor(month.name[as.numeric(as.character(bs$month))],ordered=T,levels=month.name)
bs$realm=factor(bs$realm,ordered=T,levels=c("Antarctic","Australasia","Oceania","Afrotropics","IndoMalay", "Neotropics","Palearctic","Nearctic" ))

#biomepl=melt(biomep@data,id.vars=c("id","code","biome","realm"))
#colnames(biomepl)[grep("variable",colnames(biomepl))]="month"
#biomepl$value[biomepl$value<0]=NA

p1=useOuterStrips(xyplot(mean~monthname|realm+biome,data=bs,
                 panel=function(x,y,subscripts = subscripts){
                td=bs[subscripts,]
                panel.polygon(c(td$monthname,rev(td$monthname)),c(td$meanpsd,rev(td$meanmsd)),col=grey(0.4),border=NA)
                panel.xyplot(td$monthname,td$mean,col="black",type="l",lwd=1,subscripts=subscripts)
    },scales=list(y=list(at=c(0,100),lim=c(-20,120),cex=.75,alternating=2,tck=c(0,1)),
                  x=list(at=c(1,7,12),rot=90,alternating=1)),
    ylab="Biome",xlab.top="Geographic Realm",ylab.right="MODCF (%)", xlab="Month"),
    strip=strip.custom(par.strip.text = list(cex = .7)),strip.left=strip.custom(horizontal=TRUE,par.strip.text = list(cex = .75)))
p1$par.settings$layout.widths$strip.left[1]=13
p1$par.strip.text$lines=.65
print(p1)

png("manuscript/figures/Biome_Figures.png",width=5500,height=4000,res=600,pointsize=36,bg="white")
trellis.par.set(my.theme)
print(p1)
dev.off()




