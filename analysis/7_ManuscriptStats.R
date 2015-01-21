
cld=crop(raster("data/MCD09_deriv/MCD09_meanannual_land.tif"),extent(-180,180,-60,90))
map=raster("/mnt/data/jetzlab/Data/environ/global/worldclim/bio_12.bil")
dem=raster("/mnt/data/jetzlab/Data/environ/global/worldclim/alt.bil")

dr=stack(cld,map,dem)

d=data.frame(values(dr))
#d=na.omit(d)
d2=d %.% (function(x) filter(x, complete.cases(x)))()


res=cor.test(d2$MCD09_meanannual_land,d2$alt,method="spearman",alternative="two.sided",continuity=T)

res

format(as.numeric(nrow(d)),scientific=2,digits=2)


## number of averages
r=raster("data/MCD09_deriv/meanannual.tif")
NAvalue(r)=32767

## adjust formatting to get scientific format
options(scipen =-1000)
options(digits=2)

## nubmer of non-null cells
ncell=length(na.omit(values(r)))
## number of observations (twice daily - terra and aqua)
ndays=as.integer(as.Date("30-3-2014")-as.Date("2-1-2000"))*2
## print in sci format
signif(ncell)
signif(ndays,3)
exp(log(ncell)+log(ndays))
