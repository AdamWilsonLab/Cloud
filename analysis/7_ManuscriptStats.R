
cld=crop(raster("data/MCD09_deriv/MCD09_meanannual_land.tif"),extent(-180,180,-60,90))
map=raster("/mnt/data/jetzlab/Data/environ/global/worldclim/bio_12.bil")
dem=raster("/mnt/data/jetzlab/Data/environ/global/worldclim/alt.bil")

dr=stack(cld,map,dem)

d=data.frame(values(dr))
d=na.omit(d)


res=cor.test(d$MCD09_meanannual_land,d$alt,method="spearman",alternative="two.sided",continuity=T)

res

format(as.numeric(nrow(d)),scientific=2,digits=1)
