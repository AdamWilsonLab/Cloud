### Script to download and process the NDP-026D station cloud dataset
### to validate MODIS cloud frequencies

source("setup.R")

## Data available here http://cdiac.ornl.gov/epubs/ndp/ndp026d/ndp026d.html

download=F  #download data?
## Get station locations
if(download)   system("wget -N -nd http://cdiac.ornl.gov/ftp/ndp026d/cat01/01_STID -P data/validation/NDP026D/")

st=read.table("data/validation/NDP026D/01_STID",skip=1)
colnames(st)=c("StaID","LAT","LON","ELEV","ny1","fy1","ly1","ny7","fy7","ly7","SDC","b5c")
st$lat=st$LAT/100
st$lon=st$LON/100
st$lon[st$lon>180]=st$lon[st$lon>180]-360
st=st[,c("StaID","ELEV","lat","lon")]
colnames(st)=c("id","elev","lat","lon")
write.csv(st,"stations.csv",row.names=F)
coordinates(st)=c("lon","lat")
projection(st)="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
st@data[,c("lon","lat")]=coordinates(st)

## download data
if(download){
    system("wget -N -nd ftp://cdiac.ornl.gov/pub/ndp026d/cat67_78/* -A '.tc.Z' -P data/validation/NDP026D/")
    system("gunzip data/*.Z")
}

## define FWF widths
f162=c(5,5,4,7,7,7,4) #format 162
c162=c("StaID","YR","Nobs","Amt","Fq","AWP","NC")

## use monthly timeseries
cld=do.call(rbind.data.frame,mclapply(sprintf("%02d",1:12),function(m) {
  d=read.fwf(list.files("data/validation/NDP026D/",pattern=paste("MNYDC.",m,".tc$",sep=""),full=T),skip=1,widths=f162)
  colnames(d)=c162
  d$month=as.numeric(m)
  print(m)
  return(d)}
  ))

## add lat/lon
cld[,c("lat","lon")]=coordinates(st)[match(cld$StaID,st$id),]

## drop missing values
cld=cld[,!grepl("Fq|AWP|NC",colnames(cld))]
cld$Amt[cld$Amt<0]=NA
cld$Amt=cld$Amt/100
cld=cld[!is.na(cld$Amt),]

## table of stations with > 20 observations per month
dcast(cld,StaID~YR,value.var="Nobs")
mtab=ddply(cld,c('StaID','month'),function(df){ data.frame(count=sum(df$Nobs>20,na.rm=T))})
#mtab2=mtab[table(mtab$count>10)]
stem(mtab$count)

## calculate means and sds for full record (1970-2009)
stem(cld$Nobs)
Nobsthresh=20 #minimum number of observations to include 

cldm=do.call(rbind.data.frame,by(cld,list(month=as.factor(cld$month),StaID=as.factor(cld$StaID)),function(x){
  data.frame(
      month=x$month[1],
      StaID=x$StaID[1],
      cld_all=mean(x$Amt[x$Nobs>=Nobsthresh],na.rm=T),  # full record
      cldsd_all=sd(x$Amt[x$Nobs>=Nobsthresh],na.rm=T),
      cldn_all=length(x$Amt[x$Nobs>=Nobsthresh]),
      cld=mean(x$Amt[x$YR>=2000&x$Nobs>=Nobsthresh],na.rm=T), #only MODIS epoch
      cldsd=sd(x$Amt[x$YR>=2000&x$Nobs>=Nobsthresh],na.rm=T),
      cldn=length(x$Amt[x$YR>=2000&x$Nobs>=Nobsthresh]))}))

    cldm[,c("lat","lon")]=coordinates(st)[match(cldm$StaID,st$id),c("lat","lon")]



## add the EarthEnvCloud data to cld
mod09_mean=stack(list.files("data/MCD09/",pattern="MCD09_mean_[0-9]*[.]tif",full=T))
NAvalue(mod09_mean)=255
names(mod09_mean)=month.name

mod09_sd=stack(list.files("data/MCD09/",pattern="MCD09_sd_[0-9]*[.]tif",full=T))
NAvalue(mod09_sd)=255
names(mod09_sd)=month.name


## overlay the data with 32km diameter (16km radius) buffer
## buffer size from Dybbroe, et al. (2005) doi:10.1175/JAM-2189.1.
buf16=16000
buf5=5000
bins=cut(st$lat,10)
rerun=F
if(rerun&file.exists("data/validation/valid.csv")) file.remove("data/validation/valid.csv")

beginCluster(12)

lapply(levels(bins),function(lb) {
  l=which(bins==lb)
  ## mean
  td1=raster::extract(mod09_mean,st[l,],buffer=buf5,fun=mean,na.rm=T,df=T)
  td1$id=st$id[l]
  td1$type="MCD09_meanb5"
  ## mean buffered
  td2=raster::extract(mod09_mean,st[l,],buffer=buf16,fun=mean,na.rm=T,df=T)
  td2$id=st$id[l]
  td2$type="MCD09_meanb16"
  ## std
  td3=raster::extract(mod09_sd,st[l,],buffer=buf5,fun=mean,na.rm=T,df=T)
  td3$id=st$id[l]
  td3$type="MCD09_sdb5"
  ## std buffered
  td4=raster::extract(mod09_sd,st[l,],buffer=buf16,fun=mean,na.rm=T,df=T)
  td4$id=st$id[l]
  td4$type="MCD09_sdb16"
  print(lb)#as.vector(c(l,td[,1:4])))
  write.table(rbind(td1,td2,td3,td4),"valid.csv",append=T,col.names=F,quote=F,sep=",",row.names=F)
  return(lb)
})

endCluster()

## read it back in
mod09st=read.csv("data/validation/valid.csv",header=F)[,-c(1)]
colnames(mod09st)=c(names(mod09_mean),"id","type")
mod09stl=melt(mod09st,id.vars=c("id","type"))
colnames(mod09stl)[grep("variable",colnames(mod09stl))]="month"
mod09stl$value[mod09stl$value<0]=NA
mod09stl=dcast(mod09stl,id+month~type,value="value")

## add it to cld
cldm$monthname=month.name[cldm$month]
cldm$MCD09_meanb16=mod09stl$MCD09_meanb16[match(paste(cldm$StaID,cldm$monthname),paste(mod09stl$id,mod09stl$month))]
cldm$MCD09_meanb5=mod09stl$MCD09_meanb5[match(paste(cldm$StaID,cldm$monthname),paste(mod09stl$id,mod09stl$month))]
cldm$MCD09_sdb16=mod09stl$MCD09_sdb16[match(paste(cldm$StaID,cldm$monthname),paste(mod09stl$id,mod09stl$month))]


## LULC
#system(paste("gdalwarp -r near -co \"COMPRESS=LZW\" -tr ",paste(res(mod09),collapse=" ",sep=""),
#             "-tap -multi -t_srs \"",   projection(mod09),"\" /mnt/data/jetzlab/Data/environ/global/landcover/MODIS/MCD12Q1_IGBP_2005_v51.tif ../modis/mod12/MCD12Q1_IGBP_2005_v51.tif"))
lulc=raster("/mnt/data/personal/adamw/projects/interp/data/modis/mod12/MCD12Q1_IGBP_2005_v51.tif")
require(plotKML); data(worldgrids_pal)  #load IGBP palette
IGBP=data.frame(ID=0:16,col=worldgrids_pal$IGBP[-c(18,19)],stringsAsFactors=F)
IGBP$class=rownames(IGBP);rownames(IGBP)=1:nrow(IGBP)
levels(lulc)=list(IGBP)
## function to get modal lulc value
Mode <- function(x) {
      ux <- na.omit(unique(x))
        ux[which.max(tabulate(match(x, ux)))]
      }
lulcst=raster::extract(lulc,st,fun=Mode,buffer=buf16,df=T)
colnames(lulcst)=c("id","lulc")
lulcst$StaID=st$id
## add it to cld
cldm$lulc=lulcst$lulc[match(cldm$StaID,lulcst$StaID)]
cldm$lulcc=IGBP$class[match(cldm$lulc,IGBP$ID)]


### Add biome data
biome=readOGR("data/src/teow/","biomes")
projection(biome)=projection(st)
#st$biome=over(st,biome,returnList=F)$BiomeID
dists=apply(gDistance(st,biome,byid=T),2,which.min)
st$biomec=biome$code[dists]
st$realm=biome$realm[dists]
st$biome=biome$biome[dists]

cldm$biomec=st$biomec[match(cldm$StaID,st$id)]
cldm$realm=st$relam[match(cldm$StaID,st$id)]
cldm$biome=st$biome[match(cldm$StaID,st$id)]

cldm$seas=ifelse(cldm$month%in%c(12,1,2),"DJF",
                 ifelse(cldm$month%in%3:4,"MAM",
                        ifelse(cldm$month%in%5:8,"JJA",
                               ifelse(cldm$month%in%9:11,"SON",NA))))

## identify which stations have data in which era (pre or post-MODIS)
st_era=dcast(cldml[cldml$variable%in%c("cld","cld_all"),],
             StaID~variable,value.var="value",
             fun.aggregate=function(x) sum(!is.na(x),na.rm=T))
st_era$era=ifelse(st_era$cld_all>0&st_era$cld>0,"1970-2009",
                  ifelse(st_era$cld_all>0&st_era$cld==0,"1970-2000",
                         ifelse(st_era$cld_all==0&st_era$cld>0,"2000-2009","No Data")))
st$era=st_era$era[match(st$id,st_era$StaID)]
st$era[is.na(st$era)]="No Data"

## write out the tables
write.csv(cld,file="data/validation/cld.csv",row.names=F)
write.csv(cldm,file="data/validation/cldm.csv",row.names=F)
writeOGR(st,dsn="data/validation/",layer="stations",driver="ESRI Shapefile",overwrite_layer=T)
#########################################################################
cldm=read.csv("data/validation/cldm.csv")

### Extract regional transects
cids=c(10866,10980,11130,11135,11138,11146,11210,11212,13014)

cids=c(10858,10865,10866,10980,11130,11135,11138,11204,11207,16040) #11138,11146,11210,11212,13014)
ts=cldm[cldm$StaID%in%cids,]
ts$trans=match(ts$StaID,cids)
ts$cld_allpsd=ts$cld_all+ts$cldsd_all
ts$cld_allmsd=ts$cld_all-ts$cldsd_all
ts$mod09msd=ts$mod09-ts$mod09sd
ts$mod09psd=ts$mod09+ts$mod09sd
ts$mod09msd=ts$mod09-ts$mod09sd
ts$mod09psd=ts$mod09+ts$mod09sd

tsl=melt(ts,id.vars=c("month","trans"),measure.vars=c("cld_all","cld_allmsd","cld_allpsd","MCD09_meanb5","MCD09_meanb16","MCD09_sdb16"))

#as("SpatialLinesDataFrame")
xyplot(value~trans|month,groups=variable,data=tsl,type="l",auto.key=T,ylim=c(0,100))

xyplot(cld_all~mod09|month,data=ts,type="p",auto.key=T)

