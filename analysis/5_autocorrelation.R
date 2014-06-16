source("analysis/setup.R")

#tropics=extent(c(20,25,0,5))

## mask cloud values where MAP is missing
#cf_mean=mask(cf_mean,map)

## define tropics and crop to tropics
tropics=extent(c(-180,180,-23.4378,23.4378))
global=extent(c(-180,180,-60,60))


prods=list(
  mac=raster("data/MCD09_deriv/MCD09_meanannual.tif"),
  map=raster("/mnt/data/jetzlab/Data/environ/global/worldclim/bio_12.bil"),
  dem=raster("/mnt/data/jetzlab/Data/environ/global/worldclim/alt.bil"),
  patmos=raster("data/src/gewex/CA_PATMOSX_NOAA.nc",varname="a_CA"))


region=regs[["Venezuela2"]]
regionname="Venezuela"



## loop through products and write out autocorrelation data
foreach(i=1:4) %dopar% {

tprod=names(prods[i])
  
## set file names
treg=paste0("data/autocorr/data_",tprod,"_",regionname,".tif")
tac=paste0("data/autocorr/ac_",tprod,"_",regionname,".tif")
tdist=paste0("data/autocorr/dist_",tprod,"_",regionname,".tif")

## create subset
#if(!file.exists(treg)) 
reg=crop(prods[[i]],region,filename=treg,overwrite=T,dataType='INT1S',NAflag=-128)

## mask ocean
#tcld=raster(mac_tropic)
#tmap=raster(map_tropic)
#tcld[is.na(tmap)]=NA
#tpatmos=raster(mac_tropic_patmos)

## run the autocorrelation function and write out the output raster (and time it!)
# tropics
system.time(ac<<-acorr(reg,filename=tac,gain=100,overwrite=T,dataType='INT1S',NAflag=-128))
dist=acorr_dist(reg)

#system(paste("gdal_proximity.py ",fac_trop_patmos," ",patmos_trop_dist," -co COMPRESS=LZW -co PREDICTOR=2 -ot Int16",
#             " -values 1000 -distunits GEO -nodata -32768"))


#########################################
## build the table of values to construct the correlograms
ac=raster(tac)

## summarize into tables
ftd=rbind.data.frame(
  data.frame(values=values(ac),dist=values(dist),type=tprod,region=regionname)
)
ftd$dist=round(ftd$dist)
ftd <- filter(ftd, dist <= 10000)

ftd2 <- group_by(ftd, dist,type,region)
ftd2 <- summarise(ftd2,
                  min = min(values, na.rm = TRUE),
                  max = max(values, na.rm = TRUE),
                  sd = sd(values, na.rm = TRUE),
                  mean = mean(values, na.rm = TRUE)
)

write.csv(ftd2,paste0("data/autocorr/table_",tprod,"_",regionname,".csv"),row.names=F)
print(paste("Finished ",tprod," for ",regionname))

}  ## end loop over products


## compile all regions and products into a single table

ftd3=do.call(rbind.data.frame,lapply(list.files("data/autocorr/",pattern="table",full=T),function(f) read.csv(f)))

ftdl=melt(ftd3,id.vars=c("dist","type","region"))
ftdl[,c("var","met")]=do.call(rbind,strsplit(as.character(ftdl$variable),"_"))


## plot it...
xyplot(mean~dist,data=ftd3,group=type,auto.key=F,
  panel=function(x,y,subscripts){
    td=ftd3[subscripts,]
    for(i in unique(td$type)){
    td2=td[td$type==i,]
    #  td$dist=log(td$dist+1)
    panel.xyplot(td2$dist,td2$mean,pch=16)  
    panel.segments(td2$dist,td2$min,td2$dist,td2$max,lwd=.5)
  }},
subscripts=T,xlim=c(-5,250),ylim=c(-30,104),scales=list(x=list(log=F)))


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

