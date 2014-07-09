source("analysis/setup.R")

#tropics=extent(c(20,25,0,5))

#library(devtools) 
#install_github("adammwilson/rasterAutocorr")
#detach("package:rasterAutocorr", unload=TRUE)
#library(rasterAutocorr)
## mask cloud values where MAP is missing
#cf_mean=mask(cf_mean,map)

## define tropics and crop to tropics
tropics=extent(c(-180,180,-23.4378,23.4378))
global=extent(c(-180,180,-60,60))


## create a version of the cloud data with a land mask
if(!file.exists("data/MCD09_deriv/MCD09_meanannual_land.tif"))
  system(paste0("pksetmask -i data/MCD09_deriv/meanannual.tif ",
    "-m /mnt/data/jetzlab/Data/environ/global/worldclim/alt.bil --operator='<' --msknodata -1000 --nodata 65535 ",
              "-m data/MCD09_deriv/meanannual.tif --operator='>' --msknodata 10000 --nodata 65535 ",
              " -o data/MCD09_deriv/MCD09_meanannual_land.tif"))

prods=list(
  mac=raster("data/MCD09_deriv/MCD09_meanannual_land.tif"),
  map=raster("/mnt/data/jetzlab/Data/environ/global/worldclim/bio_12.bil"),
  dem=raster("/mnt/data/jetzlab/Data/environ/global/worldclim/alt.bil"),
  patmos=raster("data/src/gewex/CA_PATMOSX_NOAA.nc",varname="a_CA"))



##### define region
region=regs[["Venezuela2"]]
regionname="Venezuela"

region=tropics
regionname="tropics"

region=regs[["SouthAmerica"]]
regionname="SouthAmerica"


## read biomes
biome=readOGR("data/src/teow/","biomes")
bcode=unique(data.frame(icode=biome$icode,code=biome$code,realm=biome$realm,biome=biome$biome))

region= biome[biome$realm=="Neotropics"&biome$biome=="Tropical & Subtropical Moist Broadleaf Forests",]

#wc=crop(prods[[2]],region)
#cld=crop(prods[[1]],region)
#dem=crop(prods[[3]],region)
#plot(stack(wc,cld,dem))
i=1
## loop through products and write out autocorrelation data
foreach(i=1:length(prods)) %do% {
tprod=names(prods[i])
acorr_table(prods[[i]],region)
write.csv(ftd2,paste0("data/autocorr/table_",tprod,"_",regionname,".csv"),row.names=F)
print(paste("Finished ",tprod," for ",regionname))

}  ## end loop over products


## compile all regions and products into a single table
#ftd3=ftd2
ftd3=do.call(rbind.data.frame,lapply(list.files("data/autocorr/",pattern=paste0("table_.*_",regionname),full=T),function(f) read.csv(f)))
ftd3$dist3=ftd3$dist2
#ftd3$dist3[ftd3$dist==0]=.01
#ftdl=melt(ftd3,id.vars=c("dist","type","region"))
#ftdl[,c("var","met")]=do.call(rbind,strsplit(as.character(ftdl$variable),"_"))

## update variable names 
levels(ftd3$type)
ftd3$type2=factor(ftd3$type,labels=c("Elevation","Mean Annual Cloud Frequency (MODCF)","Mean Annual Precipitation (WorldClim)","Mean Annual Cloud Frequency (PATMOS-X)"))


panel.bands=function(x, y, upper, lower, subscripts, ..., font, fontface) {
    upper <- upper[subscripts]
    lower <- lower[subscripts]
#    panel.polygon(c(x, rev(x)), c(upper, rev(lower)), ...)
    panel.segments(x,lower,x,upper, ...)}
x_at=c(1,2,5,10,20,50,100,200,500,1000,2000)
x_labels <- formatC(x_at, digits = 0, format = "f")

## plot it...
p1=xyplot(mean~dist3|region,data=ftd3,auto.key=list(space="inside",x=.55,y=.93),groups=type,
       upper=ftd3$mean+ftd3$sd,lower=ftd3$mean-ftd3$sd,
        panel=function(x,y,...){
            panel.superpose(x, y, panel.groups = 'panel.bands',alpha=.2,...)
            panel.xyplot(x,y,pch=16,type="l",...)  
        },
       ylab="Autocorrelation",xlab="Distance (km)",
subscripts=T,xlim=c(0.95,2000),ylim=c(-.2,1.1),scales=list(x=list(log=T,at=x_at,labels=x_labels)))+layer(panel.abline(h=0,col="red"))
p1

## save it to png
png("manuscript/figures/autoCorr.png",width=3000,height=2000,res=300,pointsize=42,bg="white")
trellis.par.set(my.theme)
print(p1)
dev.off()


## get some stats
ftd3%.%filter(dist2>10)%.% group_by(type2) %.%  arrange(desc(mean)) %.%head(1)



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

