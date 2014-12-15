
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
