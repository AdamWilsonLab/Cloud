source("analysis/setup.R")

tropics=extent(c(-180,180,-25,25))


cf_mean=raster("data/MCD09_deriv/MCD09_meanannual.tif")

tcld=crop(cf_mean,tropics)

tmap=crop(raster("/mnt/data/jetzlab/Data/environ/global/worldclim/bio_12.bil"),tcld)
projection(tmap)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

system.time(a1<<-acorr(tcld,fmean="mean"))
system.time(a2<<-acorr(tmap,fmean="mean"))

## get distances for each shift to facilitate plotting of the correlogram
d1=acorr_dist(tcld)
## get directions for each shift to facilitate plotting of the correlogram
#d2=acorr_dir(tcld)


## plot the autocorrlation and distance
plot(stack(tmap,tcld),ylab="Y",xlab="X",main="Original (Simulated) Raster")
plot(stack(a1,a2),ylab="Shift in Y",xlab="Shift in X",main="Autocorrelation")
plot(d1,ylab="Shift in Y",xlab="Shift in X",main="Distance from center in units of original raster")
plot(d2,ylab="Shift in Y",xlab="Shift in X",main="Direction from center (degrees)")


ftd=data.frame(cloud=values(a1),precip=values(a2),dist=values(d1))
ftd=ftd[ftd$dist<1000,]

## bin by dist class
ftd$distkm=round(ftd$dist)

ftdl=melt(ftd,id.vars="distkm")
ftdl2=with(ftdl,by(ftdl,by=variable,function(x) c(mean(x$value),sd(x$value))))

head(ftdl2)

xyplot(value~distkm,groups=variable,data=ftdl,auto.key=T)#[sample(1:nrow(ftd),10000),])

cellStats(dist,summary)

bwplot(value~cut(dist,c(0,5,10,20,30,40,50,100,200,1000))|variable,type=c("p","smooth"),data=ftdl,
       ylab="Correlation",xlab="Distance",main="Correlogram",sub="Different pixels within a distance class correspond to shifts of different directions (north, south, etc.) from the origin",col="black",fill="grey",scales=list(x=list(rot=45)))+
  layer(panel.abline(h=0))


