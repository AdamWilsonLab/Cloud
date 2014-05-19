source("analysis/setup.R")

cf_mean=raster("/Users/adamw/Downloads/clouds/mcd09tif_g3/MCD09_01.tif")
r="Venezuela"
tcld=crop(cf_mean,regs[[r]])
#tmap=crop(raster("/mnt/data/jetzlab/Data/environ/global/worldclim/bio_12.bil"),tcld)
#projection(tmap)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#td=stack(tmap,tcld)


x=as.matrix(tcld)

## from http://www.johnloomis.org/ece563/notes/freq/autoself/autoself.htm
library(waved)

f = fft(x);
f = f*Conj(f);
f = fft(f,inverse=T);
#f = fftshift(f);
f = Re(f);
f = f/max(max(f));
f=data.frame(f=f,id=1:length(f))
#rownames(f) = 1:nrow(f)

f2=unique(f)
f2=data.frame(id=1:length(f),r=f)
str(f)

# https://stat.ethz.ch/pipermail/r-help/2008-May/161083.html
x2=Re(fft(fft(x)* fft(x), inverse=TRUE))
image(x2)
hist(x2)
str(x2)

library(fftw)
xcorr=function(x){
  fftx=fftw2d(x)
  fftx2 = Re(fftw2d(fftx  * Conj(fftx),inverse=1))
  #im=Im(fftx2)
  #fftx2= Re(fftshift(fftx))
  return(ret)
}
str(fftx2)
hist(fftx2)
image(fftx2)

plot(fftx2[1,],type="l")

plot(f~id,data=f[sample(10000),])

n=nrow(x)*ncol(x)
x_pad = x#[x zeros(size(x))];
X     = fft(x_pad);
X_psd = abs(X)^2;
r_xx = fft(X_psd,inverse=T);
nfft = 2^nextn(2*n-1);
r = fft( fft(r_xx,nfft)* Conj(fft(r_xx,nfft)),inverse=TRUE );
# rearrange and keep values corresponding to lags: -(len-1):+(len-1)
r = [r(end-len+2:end) ; r(1:len)];



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
