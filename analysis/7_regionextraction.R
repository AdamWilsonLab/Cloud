
## for mean
mod09c=stack(list.files("data/MCD09/",pattern="MCD09_mean_[0-9]*[.]tif",full=T))  
names(mod09c)=paste0("cld_mean_",sprintf("%02d",1:12))
NAvalue(mod09c)=255
## for interannual SD:
mod09_sd=stack(list.files("data/MCD09/",pattern="MCD09_sd_[0-9]*[.]tif",full=T))       
names(mod09_sd)=paste0("cld_sd_",sprintf("%02d",1:12))
NAvalue(mod09_sd)=255

## intraannual and interannual cloud variability 
intra=raster("data/MCD09_deriv/intra.tif")
names(intra)="cld_intra"
inter=raster("data/MCD09_deriv/inter.tif")
names(inter)="cld_inter"
NAvalue(inter)=255

seas=stack("data/MCD09_deriv/seas_conc.tif")
names(seas)=c("cld_seasconc","cld_seastheta")
NAvalue(seas)=65535                                                                                                          

rdata=stack(mod09c,mod09_sd,inter,intra,seas)
names(rdata)

## read in points 
d=readOGR("data/tmp/","BK_clim_rand")

d2=raster::extract(rdata,d,sp=T)
d2@data[,c("longitude","latitude")]=coordinates(d2)

write.csv(d2@data,"data/tmp/BK_clim_cloud_rand.csv")


c(xyplot(cld_mean_06~MMP_06,data=d2@data,
         xlab="MMP_06                                                                                 Elevation"),
  xyplot(cld_mean_06~elevatn,data=d2@data),x.same=F)
