source("analysis/setup.R")

#cf_mean=raster("/Users/adamw/Downloads/clouds/mcd09tif_g3/MCD09_01.tif")
cf_mean=raster("data/MCD09_deriv/MCD09_meanannual.tif")
r="Venezuela"
tcld=cf_mean#
tcld=crop(cf_mean,regs[[r]])
#tmap=crop(raster("/mnt/data/jetzlab/Data/environ/global/worldclim/bio_12.bil"),tcld)
#projection(tmap)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#td=stack(tmap,tcld)


x=as.matrix(tcld)

## create padded version
xpad=as.matrix(extend(tcld, extent(c(-39,90,-20,30)), value=0))

library(fields)

td2=vgram.matrix(x, R=50, dx = 2,dy = 2)

system.time(td4<<-vgram.matrix(x, R=200, dx = 4,dy = 4))

save(td,file="data/out/variogram.Rdata")

####
load("data/out/variogram.Rdata")

plot(td4$d, td4$vgram, xlab="Distance")
plot(td4$d.full, td4$vgram.full, xlab="Distance")

plot(td)

head(td)



#### Use xcorr2 to get spatial autocorrelation of complete image
## http://octave.sourceforge.net/signal/function/xcorr2.html

library(RcppOctave)
o_eval("pkg load signal")

#system.time(x2<<-.O$xcorr2(x,"coeff"))
x2=.O$xcorr2(x,"coeff")

save(x2,file="data/out/meanannualxcorr2.Rdata")


load("data/out/meanannualxcorr2.Rdata")
x2=t(x2)

## get distances from center point
center=c(round(nrow(x2)/2),round(ncol(x2)/2))
dist=tcld
dist[,]=NA
dist[center[1],center[2]]=1
dist=distance(dist)

levelplot(tcld)
levelplot(dist)
levelplot(x2)

str(x2)

plot(x2[round(nrow(x2)/2),round(ncol(x2)/2):ncol(x2)],type="l",ylab="Autocorrelation",xlab="Distance (pixels)")





## from http://www.johnloomis.org/ece563/notes/freq/autoself/autoself.htm
#library(waved)

fftshift2=function(x){
  ## http://stackoverflow.com/questions/5735720/effcient-way-to-do-fft-shift-in-matlab-without-using-fftshift-function
   sz = ceiling(dim(x))/2
  x = x[
    c(sz[1]:(sz[1]*2), 1:(sz[1]-1)),
    c(sz[2]:(sz[2]*2), 1:(sz[2])-1)]
}
# https://stat.ethz.ch/pipermail/r-help/2008-May/161083.html

acorr=function(x){
  ## convert to matrix
  xm=as.matrix(x)
#  xm=x-mean(x)#xpad-mean(xpad)
  fftx=fft(xm)
  fftx2=Re(fft(fftx* Conj(fftx), inverse=TRUE))
  acor1=fftshift2(fftx2)
  acor2=acor1/max(acor1)
  res=x
  values(res)=as.vector(acor2)
  return(res)
}

a1=acorr(tcld)

## get distances from center point
center=c(round(nrow(a1)/2),round(ncol(a1)/2))
dist=raster(a1)
dist[,]=NA
dist[center[1],center[2]]=1
dist=distance(dist)/1000

td=data.frame(cor=values(a1),dist=values(dist))

plot(cor~dist,data=td)

image.plot(a1)

levelplot(tcld)
levelplot(dist)
levelplot(x2)


image.plot(a1)
image.plot(x2)

hist(x2)
str(x2)





str(r)

mcld=foreach(lag=c(3,11,51,101,201,301,501),.combine=rbind.data.frame) %dopar%
  data.frame(lag=lag,map=Moran(tcld,matrix(1,lag,lag)))

mlcdl=melt(mcld,id.vars="lag")
xyplot(value~lag,group=variable,data=mlcdl,auto.key=T,type="l")

tmap_m=MoranLocal(tmap,w=w)
tcld_m=MoranLocal(tcld,w=w)
rcor=stack(tcld_m,tmap_m)

library(fields)

histogram(rcor)
plot(rcor)
plot()
library(spatial)
kr <- krig(tmap)

pts=sampleRandom(tcld,sp=T,size=1000)
#pts=as(tmap,"SpatialPointsDataFrame")

library(pgirmess)
pgi.cor <- correlog(coords=coordinates(pts), z=pts@data[,1], method="Moran", nbclass=10)

str(pgi.cor)

#r=seq(0,30,by=.1); d=30; alpha = 0; se = 1
#kg=surf.gls(2,covmod=expcov,x=coordinates(pts)[,1],y=coordinates(pts)[,2],z=pts@data$bio_12,r=r, #=d, alpha = alpha, se = se)
#kg1=surf.gls(2,x=coordinates(pts)[,1],y=coordinates(pts)[,2],z=pts@data$MCD09_meanannual,r=r, d=d, #alpha = alpha, se = se)
#out<- mKrig( x=coordinates(pts),y=pts@data[,1],)

#image()
#str(r_xx)
#p1=correlogram(kg,nint=100,plotit=F)
#p2=correlogram(kg1,nint=100,plotit=F)

#ac1=acf(as.matrix(tmap),lag.max=10,plot=F)

#plot()

library(spdep)
# 'nb' - neighbourhood of each cell
r.nb <- dnearneigh(coordinates(pts), d1=0, d2=10)
# 'nb' - an alternative way to specify the neighbourhood
#r.nb <- cell2nb(nrow=nrow(tcld), ncol=ncol(tcld), type="queen")
sp.cor <- sp.correlogram(r.nb, pts@data[,1], order=5,
                         method="I", randomisation=T)
plot(sp.cor)

levelplot(stack(tmap_m,tcld_m))
