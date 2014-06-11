source("analysis/setup.R")

#tropics=extent(c(20,25,0,5))

cf_mean=raster("data/MCD09_deriv/MCD09_meanannual.tif")
map=raster("/mnt/data/jetzlab/Data/environ/global/worldclim/bio_12.bil")

## mask cloud values where MAP is missing
#cf_mean=mask(cf_mean,map)

## define tropics and crop to tropics
tropics=extent(c(-180,180,-23.4378,23.4378))
global=extent(c(-180,180,-60,60))

## set file names
map_global="data/autocorr/global_map.tif"
mac_global="data/autocorr/global_mac.tif"
fac_global_mac="data/autocorr/ac_global_mac.tif"
fac_global_map="data/autocorr/ac_global_map.tif"
f_global_dist="data/autocorr/ac_global_dist.tif"

map_tropic="data/autocorr/tropics_map.tif"
mac_tropic="data/autocorr/tropics_mac.tif"
fac_trop_mac="data/autocorr/ac_tropics_mac.tif"
fac_trop_map="data/autocorr/ac_tropics_map.tif"
f_trop_dist="data/autocorr/ac_tropics_dist.tif"

## create subsets
if(!file.exists(mac_global)) gcld=crop(cf_mean,global,filename=mac_global,overwrite=T,dataType='INT1S',NAflag=-128)
if(!file.exists(map_global)) crop(map,gcld,filename=map_global,overwrite=T,dataType='INT1S',NAflag=-128)

if(!file.exists(mac_tropic)) tcld=crop(cf_mean,tropics,filename=mac_tropic,overwrite=T,dataType='INT1S',NAflag=-128)
if(!file.exists(map_tropic)) crop(map,tcld,filename=map_tropic,overwrite=T,dataType='INT1S',NAflag=-128)

## read them in
gcld=raster(mac_global)
gmap=raster(map_global)

tcld=raster(mac_tropic)
tmap=raster(map_tropic)

projection(tmap)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
projection(gmap)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"


## run the autocorrelation function and write out the output raster (and time it!)
# tropics
system.time(acorr(tcld,filename=fac_tropics_mac,gain=100,overwrite=T,dataType='INT1S',NAflag=-128))
system.time(acorr(tmap,filename=fac_tropics_map,gain=100,overwrite=T,dataType='INT1S',NAflag=-128))

## global
xm=as.matrix(tcld)
tmean=mean(xm,na.rm=T)
xm=xm-tmean
## fill missing values with mean
xm[is.na(xm)]=tmean
## should we pad the image here?
## create padded version
#xpad=as.matrix(extend(tcld, extent(c(-39,90,-20,30)), value=0))
## take the fft of the matrix

object.size(x)
Rprof(tf <- "rprof.log", memory.profiling=TRUE)
fftx=fft(xm)
Rprof(NULL)
summaryRprof(tf)
system.time(fftx<<-fft(x))

system.time(g1<<-acorr(gcld,filename=fac_global_mac,gain=100,overwrite=T,dataType='INT1S',NAflag=-128))

system.time(g2<<-acorr(gmap,file=fac_global_map,gain=100,overwrite=T,dataType='INT1S',NAflag=-128))

## get distances for each shift to facilitate plotting of the correlogram
system(paste("gdal_proximity.py ",fac_global_mac," ",f_global_dist," -co COMPRESS=LZW -co PREDICTOR=2 -ot Int16",
             " -values 1000 -distunits GEO -nodata -32768"))
system(paste("gdal_proximity.py ",fac_trop_mac," ",f_trop_dist," -co COMPRESS=LZW -co PREDICTOR=2 -ot Int16",
             " -values 1000 -distunits GEO -nodata -32768"))


#########################################
## build the table of values to construct the correlograms
ac_tropics_map=raster(fac_tropics_map)
ac_tropics_mac=raster(fac_tropics_mac)
tropics_dist=raster(f_trop_dist)

ac_global_map=raster(fac_global_map)
ac_global_mac=raster(fac_global_mac)
global_dist=raster(f_global_dist)

## summarize into tables
ftd=data.frame(mac=values(ac_tropics_mac),map=values(ac_tropics_map),dist=values(tropics_dist),type="Tropics")

ftd <- filter(ftd, dist <= 5000)

ftd2 <- group_by(ftd, dist)
ftd2 <- summarise(ftd2,
                 mac_min = min(mac, na.rm = TRUE),
                  mac_max = max(mac, na.rm = TRUE),
                  mac_sd = sd(mac, na.rm = TRUE),
                  mac_mean = mean(mac, na.rm = TRUE),
                  map_min = min(map, na.rm = TRUE),
                  map_max = max(map, na.rm = TRUE),
                  map_sd = sd(map, na.rm = TRUE),
                  map_mean = mean(map, na.rm = TRUE),
                  type="tropics")

ftdl=melt(ftd2,id.vars=c("dist"))
ftdl[,c("var","met")]=do.call(rbind,strsplit(as.character(ftdl$variable),"_"))


xyplot(mac_mean~dist,data=ftd2,panel=function(x,y,subscripts){
  td=ftd2[subscripts,]
  panel.segments(td$dist,td$mac_min,td$dist,td$mac_max,lwd=.5)
  panel.xyplot(td$dist,td$mac_mean,pch=16)  
  panel.segments(td$dist,td$map_min,td$dist,td$map_max,lwd=.5,col="green")
  panel.xyplot(td$dist,td$map_mean,pch=16)  
  },xlim=c(-10,5000),ylim=c(0,104))


## plot the autocorrlation and distance
plot(stack(tmap,tcld),ylab="Y",xlab="X",main="Original Raster")
levelplot(stack(ac_tropics_mac,ac_tropics_map),ylab="Shift in Y",xlab="Shift in X",main="Autocorrelation",xlim=c(-4000,4000))
plot(d1,ylab="Shift in Y",xlab="Shift in X",main="Distance from center in units of original raster")
plot(d2,ylab="Shift in Y",xlab="Shift in X",main="Direction from center (degrees)")




cellStats(dist,summary)

bwplot(value~cut(dist,c(0,5,10,20,30,40,50,100,200,1000))|variable,type=c("p","smooth"),data=ftdl,
       ylab="Correlation",xlab="Distance",main="Correlogram",sub="Different pixels within a distance class correspond to shifts of different directions (north, south, etc.) from the origin",col="black",fill="grey",scales=list(x=list(rot=45)))+
  layer(panel.abline(h=0))

###########################
mcld=foreach(lag=c(3,11,51,101,201,301,501),.combine=rbind.data.frame) %dopar%
  data.frame(lag=lag,map=Moran(tcld,matrix(1,lag,lag)))

mlcdl=melt(mcld,id.vars="lag")
xyplot(value~lag,group=variable,data=mlcdl,auto.key=T,type="l")

w=matrix(1,nrow=3,ncol=3)
tmap_m=MoranLocal(tmap,w=w)
tcld_m=MoranLocal(tcld,w=w)
rcor=stack(tcld_m,tmap_m)

