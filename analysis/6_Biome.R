### Produce summary of cloud frequency by biome

source("analysis/setup.R")


## Create a simplified version of the TEOW biome dataset
if(!file.exists("data/src/teow/biomes.shp")){
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
  writeOGR(biome,"data/src/teow","biomes",driver="ESRI Shapefile",overwrite=T)
}

## rasterize the biome dataset
## create a copy with a numeric biome code
biome=readOGR("data/src/teow/","biomes")
bcode=unique(data.frame(icode=biome$icode,code=biome$code,realm=biome$realm,biome=biome$biome,area=biome$areakm))

system(paste("gdal_rasterize -a icode -init 0 -l biomes -ot Byte -te -180 -90 180 90 -tr 0.008333333333333 -0.008333333333333",
              " -co COMPRESS=LZW -co ZLEVEL=9 -co PREDICTOR=2 ",
              " data/src/teow/biomes.shp  data/out/teow/teow.tif"))

#########################################
### Summarize the cloud data by biome

## create a list of products to summarize
bprods=c("data/MCD09_deriv/inter.tif",
         "data/MCD09_deriv/intra.tif",
         "data/MCD09_deriv/seasconc.tif",
         "data/MCD09_deriv/seastheta.tif",
         "data/MCD09_deriv/meanannual.tif",
    paste("data/MCD09/MCD09_mean_",sprintf("%02d",1:12),".tif",sep=""))

### loop over products and summarize by biome
foreach(m=bprods)%dopar%{
  ## set file names
  tbiome=paste0("data/tmp/teow_",basename(m))
  tcloudbiome=paste0("data/out/biomesummaries/teow_",sub(".tif",".txt",basename(m)))

  ## mask biome raster using missing data in cloud dataset
  nas=sub("^.*=","",system(paste0("gdalinfo ",m," | grep NoData"),intern=T))  #get NA for image
  system(paste("pksetmask -i data/out/teow.tif -m ",m," -ot UInt16 ",
               "--operator='=' --msknodata ",nas," --nodata 0  -co COMPRESS=LZW -co PREDICTOR=2 -o ",tbiome))
  ## calculate biome-level summary metrics
  system(paste("oft-stat -i ",m," -o ",tcloudbiome," -um ",tbiome," -mm"))
  ## clean up
  file.remove(tbiome)
}    

bs=do.call(rbind.data.frame,lapply(bprods,function(m){
  tcloudbiome=paste0("data/out/biomesummaries/teow_",sub(".tif",".txt",basename(m)))
  print(tcloudbiome)
  td=read.table(tcloudbiome)
  colnames(td)=c("icode","n","min","max","mean","sd")
  td$product=sub(".tif","",basename(m))
  td$meanpsd=td$mean+td$sd
  td$meanmsd=td$mean-td$sd
  td=merge(td,bcode,by="icode")
#  file.remove(tcloudbiome)
  return(td)
  }))
write.csv(bs,file="data/out/biomesummary.csv",row.names=F)

###################################################################
### summary by biome
bs=read.csv(file="data/out/biomesummary.csv")
bs$realm=factor(bs$realm,ordered=T,levels=c("Antarctic","Australasia","Oceania","Afrotropics","IndoMalay", "Neotropics","Palearctic","Nearctic" ))


## filter to get monthly timeseries for plotting
bsm=bs[grep("MCD09",bs$product),]
bsm$monthname=factor(month.name[as.numeric(as.character(sub("MCD09_mean_","",bsm$product)))],ordered=T,levels=month.name)

#biomepl=melt(biomep@data,id.vars=c("id","code","biome","realm"))
#colnames(biomepl)[grep("variable",colnames(biomepl))]="month"
#biomepl$value[biomepl$value<0]=NA

### plot by month
p1=useOuterStrips(xyplot(I(mean/100)~monthname|realm+biome,data=bsm,
                 panel=function(x,y,subscripts = subscripts){
                td=bsm[subscripts,]
                panel.polygon(c(td$monthname,rev(td$monthname)),c(td$meanpsd/100,rev(td$meanmsd/100)),col=grey(0.4),border=NA)
                panel.xyplot(td$monthname,td$mean/100,col="black",type="l",lwd=1,subscripts=subscripts)
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

####
## Export summary table
#bs$meansd=paste(round(bs$mean/100,1)," (",round(bs$sd/100,1),")",sep="")
#bst=dcast(bs[bs$product%in%c("meanannual","inter","intra"),],realm+biome~product,value.var="meansd")
#colnames(bst)=c("Realm","Biome","Interannual","Intraannual","MeanAnnual")


## Summary stats for paper
## overall max
levels(bs$product)

bs%.%filter(product=="meanannual")%.% arrange(desc(mean)) %.%head(5)

bs%.%filter(product=="intra")%.% arrange(desc(mean)) %.%head(5)
bs%.%filter(product=="intra")%.% arrange(mean) %.%head(5)

bs%.%filter(product=="inter")%.% arrange(desc(mean)) %.%head(5)
bs%.%filter(product=="inter")%.% arrange(mean) %.%head(5)

bs%.%filter(product=="seasconc")%.% arrange(desc(mean)) %.%head(5)
bs%.%filter(product=="seastheta")%.% arrange(desc(mean)) %.%head(5)


bs[which.max(bs$mean),]

bs[which.max(bs$min),]
bs[which.max(bs$seasintra),]


bs[which.min(bs$mean),]



