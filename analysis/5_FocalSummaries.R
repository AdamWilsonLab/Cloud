## Perform focal analysis on Cloud outputs in parallel

# Load the libraries and set working directory
source("analysis/setup.R")


### Build the tiles to process
jobs=tilebuilder(xmin=-180,xmax=180,ymin=-90,ymax=90,size=20,overlap=2)
#jobs=tilebuilder(xmin=-74,xmax=-43,ymin=-30,ymax=14,size=10,overlap=2)

#ggplot(jobs)+
#  geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="transparent",lty="dashed",col="darkgreen")+
#  geom_rect(aes(xmin=xminb,xmax=xmaxb,ymin=yminb,ymax=ymaxb),fill="transparent",col="black")+
#  geom_text(aes(x=(xminb+xmax)/2,y=(yminb+ymax)/2,label=tile),size=3)+
#  coord_equal()


## 100km spatial SD
#dmeanannual=raster("data/MCD09_deriv/meanannual.tif")
file="data/MCD09_deriv/meanannual.tif"
## get diameter of circle (in pixels) for a 1 degree circle
#fw=focalWeight(raster(file), d=0.5, type='circle')
#dim(fw)

if(!file.exists(paste0(datadir,"/mcd09focal"))) dir.create(paste0(datadir,"/mcd09focal"))

registerDoMC(20)

foreach( i=1:nrow(jobs), .options.multicore=list(preschedule=FALSE)) %dopar% {
  
  toutfile1=paste(datadir,"/mcd09focal/", sub(".tif","",basename(file)),"_",jobs$tile[i],"_region.tif",sep="")
  toutfile2=paste(datadir,"/mcd09focal/", sub(".tif","",basename(file)),"_",jobs$tile[i],"_sd.tif",sep="")
  toutfile3=paste(datadir,"/mcd09focal/", sub(".tif","",basename(file)),"_",jobs$tile[i],"_sdcrop.tif",sep="")
  
  if(file.exists(toutfile2)) {writeLines(paste(toutfile,"Exists, moving on"));return(NULL)}
  writeLines(paste("Starting: ",basename(toutfile1)," tile:",jobs$tile[i]," ( ",i," out of ",nrow(jobs),")"))
  ## crop to buffered region
  system(paste("gdalwarp -ot Int16 -srcnodata 65535 -dstnodata -32768 -te ",paste(select(jobs,xminb,yminb,xmaxb,ymaxb)[i,],collapse=" "),file,toutfile1))
  ## get focal SD
 system(paste("pkfilter -nodata -32768 -nodata 32767 -dx 121 -dy 121 -f stdev -circ TRUE -i ",toutfile1," -o ",toutfile2))
#   system(paste("pkfilter -nodata -32768 -nodata 32767 -dx 21 -dy 21 -f stdev -circ TRUE -i ",toutfile1," -o ",toutfile2))
  
  ## crop to main tile
  system(paste("gdalwarp -te ",paste(select(jobs,xmin,ymin,xmax,ymax)[i,],collapse=" "),toutfile2,toutfile3))
  ## remove temporary files
  file.remove(toutfile1,toutfile2)
  print(paste("Finished Temporary File: ",basename(toutfile1)))
}


system(paste0("gdalwarp ",datadir,"/mcd09focal/*_sdcrop.tif data/MCD09_deriv/mean_1deg_sd_uncompressed.tif"))
system(paste("gdal_translate -co COMPRESS=DEFLATE -co ZLEVEL=9 -co BIGTIFF=YES -co COMPRESS=LZW -co PREDICTOR=2",
             " data/MCD09_deriv/mean_1deg_sd_uncompressed.tif data/MCD09_deriv/mean_1deg_sd.tif"))

file.copy("data/MCD09_deriv/mean_1deg_sd.tif","data/MCD09_EarthEngineUpload/mean_1deg_sd.tif")

