###  Script to compile the monthly cloud data from earth engine into a netcdf file for further processing

library(raster)
library(doMC)
registerDoMC(12)


## start raster cluster
#beginCluster(5)

setwd("~/acrobates/adamw/projects/cloud")

datadir="/mnt/data2/projects/cloud/"


### Download files from google drive
## This only works if google-cli is installed and has already been authenticated 
download=F
if(download) system(paste("google docs get 2014043_g3_* ",datadir,"/mcd09ee",sep=""))


##  Get list of available files
version="g3"  #which version of data from EE?
df=data.frame(path=list.files(paste(datadir,"mcd09ee",sep="/"),pattern=paste(".*",version,".*.tif$",sep=""),full=T,recur=T),stringsAsFactors=F)
df[,c("month","sensor")]=do.call(rbind,strsplit(basename(df$path),"_|[.]|-"))[,c(5,4)]
df$date=as.Date(paste(2013,"_",df$month,"_15",sep=""),"%Y_%m_%d")


## use ramdisk?
tmpfs="tmp/"#tempdir()

ramdisk=F
if(ramdisk) {
    system("sudo mkdir -p /mnt/ram")
    system("sudo mount -t ramfs -o size=30g ramfs /mnt/ram")
    system("sudo chmod a+w /mnt/ram")
    tmpfs="/mnt/ram"
}

rasterOptions(tmpdir=tmpfs,overwrite=T, format="GTiff",maxmemory=1e9)
rerun=T  # set to true to recalculate all dates even if file already exists

## define month-sensors to process
jobs=unique(data.frame(month=as.numeric(df$month),sensor=df$sensor))

i=1
#jobs=jobs[jobs$sensor=="MYD09",]



### add boundaries to file list to remove problematic pixels at high latitudes
## these boundaries were later added to the earth engine script, so if it is re-run this is not necessary
#xmin,xmax,ymin,ymax
mextent=rbind.data.frame(
    "01"=c(-180,-90,180,73.5),#
    "02"=c(-180,-90,180,84),  #
    "03"=c(-180,-90,180,90),
    "04"=c(-180,-90,180,90),
    "05"=c(-180,-69,180,90),
    "06"=c(-180,-62.5,180,90),  #
    "07"=c(-180,-67,180,90),    #
    "08"=c(-180,-77,180,90),
    "09"=c(-180,-77,180,90),
    "10"=c(-180,-90,180,89), #
    "11"=c(-180,-90,180,77),
    "12"=c(-180,-90,180,69)
)    
colnames(mextent)=c("xmin","ymin","xmax","ymax")

## project to sinusoidal
proj="'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'"

mextentsin=data.frame(t(apply(mextent,1,function(x) c(project(t(x[1:2]),sub("'","",proj)),project(t(x[3:4]),sub("'","",proj))))))
colnames(mextentsin)=c("xmin","ymin","xmax","ymax")
mextentsin$xmin=project(t(c(-180,0)),sub("'","",proj))[1]
mextentsin$xmax=project(t(c(180,0)),sub("'","",proj))[1]

## Loop over data to mosaic tifs, compress, and add metadata
    foreach(i=1:nrow(jobs)) %dopar% {
        ## get month
        m=jobs$month[i]
        cm=sprintf("%02d",m)
        
        date=df$date[df$month==m][1]
        print(date)
        ## get sensor
        s=jobs$sensor[i]
        s2=sub("GA","",s)

        ## Define output and check if it already exists
        tvrt=paste(datadir,"/mcd09tif/",s2,"_",cm,"_globalsin.vrt",sep="")
        tvrt2=paste(datadir,"/mcd09tif/",s2,"_",cm,"_globalwgs84.vrt",sep="")
        ttif=paste(datadir,"/mcd09tif/",s2,"_",cm,"_mean.tif",sep="")
        ttif2=paste(datadir,"/mcd09tif/",s2,"_",cm,"_sd.tif",sep="")

        ## check if output already exists
        if(!rerun&file.exists(ttif)) return(NA)
        ## build VRT to merge tiles
        ## include subseting using mextentsin object created above to ensure cropping of problematic values
        system(paste("gdalbuildvrt -b 1 -b 2 -te ",paste(mextentsin[m,],collapse=" ")," -srcnodata -32768 -vrtnodata 32767 ",tvrt," ",paste(df$path[df$month==m&df$sensor==s],collapse=" ")))
        ## Merge to geotif in temporary directory
        ## specify sourc projection because it gets it slightly wrong by default 
        ops=paste("-multi -of vrt --config GDAL_CACHEMAX 500 -wm 500 -wo NUM_THREADS:10 -srcnodata 32767 -dstnodata 32767 -s_srs ",proj,"  -t_srs 'EPSG:4326' ",
            " -ot Int16 -r bilinear -te -180 -90 180 90 -tr 0.008333333333333 -0.008333333333333")
        ## create the warpped VRT
        system(paste("gdalwarp -overwrite ",ops," ",tvrt," ",tvrt2))

        ## Compress file and add metadata tags
        ops2=paste("-ot Int16 -co COMPRESS=LZW -co PREDICTOR=2 -stats -a_nodata 32767 ")
        meantags=c(paste("TIFFTAG_IMAGEDESCRIPTION='Monthly Cloud Frequency (x10000) for 2000-2013 extracted from C5 MODIS ",s,"GA PGE11 internal cloud mask algorithm (embedded in state_1km bit 10).",
            "The daily cloud mask time series were summarized to mean cloud frequency (CF) by calculating the proportion of cloudy days. '"),
              paste("TIFFTAG_DOCUMENTNAME='Collection 5 ",s," Mean Cloud Frequency'",sep=""),
              paste("TIFFTAG_DATETIME='2013",sprintf("%02d", m),"15'",sep=""),
              "TIFFTAG_ARTIST='Adam M. Wilson (adam.wilson@yale.edu)'")
        sdtags=c(paste("TIFFTAG_IMAGEDESCRIPTION='Standard Deviation (x10000) of Monthly Cloud Frequency for 2000-2013 extracted from C5 MODIS ",s,"GA PGE11 internal cloud mask algorithm (embedded in state_1km bit 10).",
            "The daily cloud mask time series were summarized to mean cloud frequency (CF) by calculating the proportion of cloudy days. '"),
              paste("TIFFTAG_DOCUMENTNAME='Collection 5 ",s," SD Cloud Frequency'",sep=""),
              paste("TIFFTAG_DATETIME='2013",sprintf("%02d", m),"15'",sep=""),
              "TIFFTAG_ARTIST='Adam M. Wilson (adam.wilson@yale.edu)'")

        ## run the merge, warp, and compress all in one step...
        system(paste("gdal_translate -b 1  ",ops2," ",paste("-mo ",meantags,sep="",collapse=" ")," ",tvrt2," ",ttif))
        system(paste("gdal_translate -b 2  ",ops2," ",paste("-mo ",sdtags,sep="",collapse=" ")," ",tvrt2," ",ttif2))
        writeLines(paste("Month:",m," Sensor:",s," Finished"))
    }



