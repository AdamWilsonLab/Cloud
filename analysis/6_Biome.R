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
         "data/MCD09_deriv/mean_1deg_sd.tif",
    paste("data/MCD09/MCD09_mean_",sprintf("%02d",1:12),".tif",sep=""))

### loop over products and summarize by biome
foreach(m=bprods)%do%{
  ## set file names
  tcloud=paste0("data/tmp/rescale_",sub("tif","vrt",basename(m)))
  tbiome=paste0("data/tmp/teow_",basename(m))
  tcloudbiome=paste0("data/out/biomesummaries/teow_",sub(".tif",".txt",basename(m)))
  tcloudhist=paste0("data/out/biomesummaries/teowhist_",sub(".tif",".txt",basename(m)))
  
  ## mask biome raster using missing data in cloud dataset
  nas=sub("^.*=","",system(paste0("gdalinfo ",m," | grep NoData"),intern=T))  #get NA for image
  if(grepl("1deg_sd",m)) nas=0
  system(paste("pksetmask -i data/out/teow.tif -m ",m," -ot Byte ",
               "--operator='=' --msknodata ",nas," --nodata 0  -co COMPRESS=LZW -co PREDICTOR=2 -o ",tbiome))
##select only one biome
   system(paste("pksetmask -i data/out/teow.tif -ot Byte  ",   
                " -m data/out/teow.tif --operator='=' --msknodata 18 --nodata 70 ",
                " -m data/out/teow.tif --operator='=' --msknodata 29 --nodata 70 ",
                " -m data/out/teow.tif --operator='=' --msknodata 65 --nodata 70 ",
                " -m data/out/teow.tif --operator='=' --msknodata 69 --nodata 70 ",
                " -m ",m," --operator='=' --msknodata ",nas," --nodata 0 ",
#                " -m ",m," --operator='>' --msknodata 10000 --nodata 0 ",
                " -co COMPRESS=LZW -co PREDICTOR=2 -o ",tbiome))

  
  ## calculate biome-level summary metrics
  system(paste("oft-stat -i ",m," -o ",tcloudbiome," -um ",tbiome," -mm"))
  ## calculate biome-level histograms for monthly data
  system(paste("gdal_translate  -scale 0 10000 0 100 -ot Byte -of vrt  ",m,tcloud))
  system(paste("oft-his -i ",tcloud," -o ",tcloudhist," -um ",tbiome," -hr -maxval 100 "))
               
               ## clean up
  file.remove(tbiome,tcloud)
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

## summarize histograms
fhist=list.files("data/out/biomesummaries/",pattern="teowhist_MCD09_mean_.*txt",full=T)
bsf=do.call(rbind.data.frame,lapply(fhist,
            function(tcloudbiome){
  print(tcloudbiome)
  td=read.table(tcloudbiome)
  colnames(td)=c("icode","n","band",paste("v",0:100,sep=""))
  td$product=sub(".tif","",basename(tcloudbiome))
  td=merge(td,bcode,by="icode")
  td$month=as.numeric(sub(".txt","",sub("^.*teowhist_MCD09_mean_","",tcloudbiome)))
  #  file.remove(tcloudbiome)
  return(td)
}))
write.csv(bsf,file="data/out/biomesummaryhist.csv",row.names=F)


###################################################################
### summary by biome
bs=read.csv(file="data/out/biomesummary.csv")
bs$realm=factor(bs$realm,ordered=T,levels=c("Antarctic","Australasia","Oceania","Afrotropics","IndoMalay", "Neotropics","Palearctic","Nearctic" ))

bsf=read.csv(file="data/out/biomesummaryhist.csv")
bsfl=melt(select(bsf,c(-product,-band)),id.vars=c("icode","n","code","realm","biome","month","area"),value.name="count")
bsfl$value=as.numeric(sub("v","",as.character(bsfl$variable)))
## convert area to percentage
bsfl$ncount=bsfl$count/bsfl$n


####
## Export summary table
#bs$meansd=paste(round(bs$mean/100,1)," (",round(bs$sd/100,1),")",sep="")
#bst=dcast(bs[bs$product%in%c("meanannual","inter","intra"),],realm+biome~product,value.var="meansd")
#colnames(bst)=c("Realm","Biome","Interannual","Intraannual","MeanAnnual")


## Summary stats for paper
## overall max
levels(bs$product)

bs%.%filter(product=="meanannual")%.% arrange(desc(mean)) %.%head(5)
bs%.%filter(product=="meanannual")%.% arrange(mean) %.%head(5)

bs%.%filter(product=="intra")%.% arrange(desc(mean)) %.%head(15)
bs%.%filter(product=="intra")%.% arrange(mean) %.%head(5)

bs%.%filter(product=="inter")%.% arrange(desc(mean)) %.%head(5)
bs%.%filter(product=="inter")%.% arrange(mean) %.%head(5)

bs%.%filter(product=="seasconc")%.% arrange(desc(mean)) %.%head(5)
bs%.%filter(product=="seastheta")%.% arrange(desc(mean)) %.%head(5)

## Spatial
bs%.%filter(product=="mean_1deg_sd")%.% arrange(desc(mean)) %.%head(5)
bs%.%filter(product=="mean_1deg_sd")%.% arrange(mean) %.%head(5)


bs[which.max(bs$mean),]

bs[which.max(bs$min),]
bs[which.max(bs$seasintra),]


bs[which.min(bs$mean),]


### Convert to quantiles
fquantile=function(vals,freq,quant) {
  ord <- order(vals)
  freq2=freq/sum(freq)
  cs <- cumsum(freq2[ord])
  tx=do.call(c,lapply(quant,function(tquant) vals[max(which(cs<tquant))+1] ))
  names(tx)=paste0("Q",quant)
  return(tx)
}

#fquantile(vals=bsfl$value,freq=bsfl$count,c(0.0000000000001,0.025,0.25,0.5,0.75,.975,1))
group_by(bsfl,biome)%.%summarize(area=sum(log(area)))

paste(levels(bsfl$biome),collapse="','")
biomebin=matrix(c(
  'Boreal Forests/Taiga',                                    'Other',
  'Deserts & Xeric Shrublands',                 'Deserts &\nXeric Shrublands',    
  'Flooded Grasslands & Savannas',                         'Other',
  'Lake',                                                   'Other',
  'Mangroves',                                              'Other',
  'Mediterranean Forests, Woodlands & Scrub',  'Mediterranean Forests,\nWoodlands & Scrub',
  'Montane Grasslands & Shrublands',          'Montane Grasslands &\nShrublands',
  'Rock & Ice',                                             'Other',
  'Temperate Broadleaf & Mixed Forests',      'Temperate\nForests',
  'Temperate Conifer Forests',           'Temperate\nForests',
  'Temperate Grasslands, Savannas & Shrublands', 'Temperate Grasslands,\nSavannas,\n& Shrublands',
  'Tropical & Subtropical Coniferous Forests','Tropical & Subtropical\nConiferous and\nDry Broadleaf Forests',
  'Tropical & Subtropical Dry Broadleaf Forests','Tropical & Subtropical\nConiferous and\nDry Broadleaf Forests',
  'Tropical & Subtropical Grasslands, Savannas & Shrublands', 'Tropical & Subtropical\nGrasslands, Savannas,\n& Shrublands',
  'Tropical & Subtropical Moist Broadleaf Forests', 'Tropical & Subtropical\nMoist Broadleaf\nForests',
  'Tundra',                                                   'Other'),ncol=2,byrow=T) 
colnames(biomebin)=c("old","new")

bsfl$biome2= biomebin[match(bsfl$biome,biomebin[,"old"]),"new"]


## calculate quantiles by biome
qs=group_by(bsfl,biome2,realm,month)%.% summarize(
  Q0=min(value[count>0],na.rm=T),
  Q02.5=fquantile(vals=value,freq=count,0.025),
  Q25=fquantile(vals=value,freq=count,0.25),
  Q50=fquantile(vals=value,freq=count,0.5),
  Q75=fquantile(vals=value,freq=count,0.75),
  Q97.5=fquantile(vals=value,freq=count,0.975),
  Q100=fquantile(vals=value,freq=count,1),
  areakm=sum(area[variable=="v0"])) 

pbiome=

  ggplot(qs) +
  geom_ribbon(mapping=aes(x=month, ymin=Q0,ymax=Q100),fill=grey(.7))+
  geom_ribbon(mapping=aes(x=month, ymin=Q02.5,ymax=Q97.5),fill=grey(.5))+
  geom_ribbon(mapping=aes(x=month, ymin=Q25,ymax=Q75),fill=grey(.2))+
  geom_line(mapping=aes(x=month, y=Q50),col="red")+
  geom_point(aes(x=12,y=0,size = areakm),colour="darkblue")+
  scale_size(trans="log10",guide=
               guide_legend(title="Area (km)",direction="horizontal",label.position="bottom"))+
  facet_grid(biome2~realm)+
  scale_y_continuous(breaks=seq(0, 100, 50))+
  scale_x_continuous(breaks=c(3,6,9))+
  theme_bw()+
  theme(strip.text.x = element_text(angle = 0), plot.margin=unit(c(.05,.05,.1,.05),"npc"))+
  theme(strip.text.y = element_text(angle = 0),legend.position=c(1.05,-.1))+
  ylab("Cloud Frequency (%)")+
  xlab("Month")

png(file=paste0("figure/biome_overview.png"),width=3000,height=3000,pointsize=24,res=300)
print(pbiome)
dev.off()



### Spatial plot using shapefile
bsw=dcast(bs,code~product,value.var="mean")
biome2=biome
biome2@data[,colnames(bsw)[-1]]=bsw[match(biome2$code,bsw$code),-1]

spplot(biome2,zcol="meanannual")