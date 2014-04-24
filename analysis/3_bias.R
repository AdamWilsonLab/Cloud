###  Script to compile the monthly cloud data from earth engine into a netcdf file for further processing

library(rasterVis)
library(multicore)
library(doMC)
library(foreach)
library(rgdal)
registerDoMC(12)


# final output will be written to data directory here:
setwd("/mnt/data/personal/adamw/projects/cloud")

# temporary files will be written here:
datadir="/mnt/data2/projects/cloud/"


## Specify path to VSNR souce code and add it to RcppOctave path
library(RcppOctave)
mpath="/mnt/data/personal/adamw/projects/environmental-layers/climate/research/cloud/MOD09/vsnr/"
.O$addpath(mpath)


#########################################
####  Bias correction functions

fgabor=function(d,theta,x=200,y=5){
    thetaR=(theta*pi)/180
    cds=expand.grid(x=(1:nrow(d))-round(nrow(d)/2),y=(1:ncol(d))-round(ncol(d)/2))
    sigma_x=x
    sigma_y=y
    lambda=0
    tpsi=0
    x_theta=cds[,"x"]*cos(thetaR)+cds[,"y"]*sin(thetaR);
    y_theta=-cds[,"x"]*sin(thetaR)+cds[,"y"]*cos(thetaR);
    n=max(cds)
    gb= 1/(2*pi*sigma_x *sigma_y) * exp(-.5*(x_theta^2/sigma_x^2+y_theta^2/sigma_y^2))*cos(2*pi/n*lambda*x_theta+tpsi);
    gb2=1e-2*gb/max(gb); #Normalization
    psi=d
    values(psi)=matrix(gb2,ncol=ncol(d))
    return(psi)
}
  
     
vsnr=function(d,gabor,alpha=2,p=2,epsilon=0,prec=5e-3,maxit=100,C1=1,full=F){
    ## VSNR can't run with any missing values, set them to zero here then switch them back to NA later
    d2=as.matrix(d)
    d2[is.na(d2)]=0
    ## Process with VSNR
    dt=.CallOctave("VSNR",d2,epsilon,p,as.matrix(gabor),alpha,maxit+1,prec,C1,argout ="u",verbose=F);
    ## Create spatial objects from VSNR output
    dc=d
    values(dc)=as.numeric(t(dt$u))  # Make 'corrected' version of data
    ##  Set NA values in original data back to NA 
    dc[is.na(d)]=NA
    ## remove temp files
    rm(d2)
    ## return the corrected data
    return(dc)
}

rmr=function(x){
    ## function to truly delete raster and temporary files associated with them
        if(class(x)=="RasterLayer"&grepl("^/tmp",x@file@name)&fromDisk(x)==T){
            file.remove(x@file@name,sub("grd","gri",x@file@name))
            rm(x)
    }
}

                    
######################################
## Run the correction functions
###  Subset equitorial region to correct orbital banding


### build table of tiles to process
extf=function(xmin=-180,xmax=180,ymin=-30,ymax=30,size=10,overlap=0.5){
    xmins=unique(sort(c(seq(xmin,xmax-size,by=size),seq(xmin+(overlap*size),xmax-size,by=size))))
    ymins=unique(sort(c(seq(ymin,ymax-size,by=size),seq(ymin+(overlap*size),ymax-size,by=size))))
    exts=expand.grid(xmin=xmins,ymin=ymins)
    exts$ymax=exts$ymin+size
    exts$xmax=exts$xmin+size
    exts$tile=1:nrow(exts)
    return(exts)
}


i=1
ti=7
### Build the tiles to process
exts=extf(xmin=-180,xmax=180,ymin=-30,ymax=30,size=60,overlap=0)

## add an extra tile to account for regions of reduced data availability for each sensor
modexts=cbind.data.frame(sensor="MOD09",rbind.data.frame(
                             exts,
                             c(xmin=130,ymin=-50,ymax=-25,xmax=180,tile=nrow(exts)+1),
                             c(xmin=130,ymin=-25,ymax=0,xmax=180,tile=nrow(exts)+2)))
mydexts=cbind.data.frame(sensor="MYD09",rbind.data.frame(
                             exts,
                             c(xmin=-170,ymin=-30,ymax=0,xmax=-140,tile=nrow(exts)+1),
                             c(xmin=-170,ymin=0,ymax=55,xmax=-140,tile=nrow(exts)+2)))

allexts=rbind.data.frame(modexts,mydexts)

### assemble list of files to process
df2=data.frame(path=list.files(paste(datadir,"/mcd09tif",sep=""),full=T,pattern="[0-9]*[mean|sd].*tif$"),stringsAsFactors=F)
df2[,c("sensor","month","type")]=do.call(rbind.data.frame,strsplit(basename(df2$path),"_|[.]"))[,c(1,2,3)]

## create a list of jobs to process
jobs=data.frame(allexts,month=rep(sprintf("%02d",1:12),each=nrow(allexts)))
jobs=rbind.data.frame(cbind.data.frame(type="mean",jobs),cbind.data.frame(type="sd",jobs))
jobs$path=df2$path[match(paste(jobs$sensor,jobs$month,jobs$type),paste(df2$sensor,df2$month,df2$type))]
## drop any jobs with no associated files
jobs=jobs[!is.na(jobs$path),]


## loop over sensor-months to create full grid of corrected values
foreach( i=1:nrow(jobs), .options.multicore=list(preschedule=FALSE)) %dopar% {
    file=jobs$path[i]
    toutfile=paste(datadir,"mcd09bias/", sub(".tif","",basename(file)),"_",jobs$tile[i],".tif",sep="")
    if(file.exists(toutfile)) {writeLines(paste(toutfile,"Exists, moving on"));return(NULL)}
    writeLines(paste("Starting: ",toutfile," tile:",jobs$tile[i]," ( ",i," out of ",nrow(jobs),")"))
    ## set angle of orbital artefacts to be corrected
    sensor=jobs$sensor[i]
    if(sensor=="MOD09") scanangle=-15
    if(sensor=="MYD09") scanangle=15
    ## Process the tiles
    textent=extent(jobs$xmin[i],jobs$xmax[i],jobs$ymin[i],jobs$ymax[i])
    ## extract the tile
    ## extract the data
    d=crop(raster(file),textent)
    ## acount for scale of data is 10000*CF, so convert to 0-100
    d2=d*.01
    ## skip null tiles - will only have this if tiles are quite small (<10 degrees)
    if(is.null(d2@data@values)) return(NULL)
    ## make the gabor kernel
    ## this specifies the 'shape' of the noise we are trying to remove
    psi=fgabor(d2,theta=scanangle,x=400,y=4) #3
    ## run the correction function.  
    res=vsnr(d2,gabor=psi,alpha=2,p=2,epsilon=1,prec=5e-6,maxit=50,C=1,full=F)
    res2=res*100
    ## write the file
    writeRaster(res2,file=toutfile,overwrite=T,datatype='INT2S',options=c("COMPRESS=LZW", "PREDICTOR=2"),NAvalue=32767)
    ## remove temporary files
    rmr(d);rmr(d2);rmr(psi);rmr(res);rmr(res2)
    ## remove old temporary files older than x hours
    removeTmpFiles(1)
    print(paste("Finished Temporary File: ",toutfile))
}


############################################
## now mosaic the tiles with the original image to keep only the corrected data (when available) and the uncorrected data where there is no tile.
## this relies on df2 created above
#tfs=list.files(paste(datadir,"/mcd09bias",sep=""),pattern="M.*_[0-9]._[0-9][.]tif",full=T)
#file.rename(tfs,paste(substr(tfs,1,46),"mean_",substr(tfs,47,52),sep=""))

## define color scale for mean and sd
mkct=function(palette,vals,bit)
    data.frame(id=0:(2^bit-1),t(col2rgb(c(palette(length(vals)+1),rep("#00000000",(2^bit)-length(vals)-1)),alpha=T)))

colR=colorRampPalette(c("#08306b","#0d57a1","#2878b8","#4997c9","#72b2d7","#a2cbe2","#c7dcef","#deebf7","#f7fbff"))
meancols=mkct(colR,vals=0:10000,bit=16)
write.table(meancols,file="data/colors16_mean.txt",quote=F,row.names=F,col.names=F)

colR2=colorRampPalette(c("#0000FF","#00FF80","#FF0080"))
sdcols=mkct(colR2,vals=0:10000,bit=16)
write.table(sdcols,file="data/colors16_sd.txt",quote=F,row.names=F,col.names=F)


foreach( i=1:nrow(df2)) %dopar% {
    ifile=df2$path[i]
    itype=df2$type[i]
    outfile=paste(datadir,"/mcd09ctif/",sub(".tif",".vrt",basename(ifile)),sep="")
    outfile2=paste(datadir,"/mcd09ctif/unmasked_",basename(ifile),sep="")
    outfile3=paste(datadir,"/mcd09ctif/",basename(ifile),sep="")
#    if(file.exists(outfile3)) {print(paste(outfile," exists, moving on...")); return(NULL)}
    ## mosaic the tiles with the original data (keeping the new data when available)
    tfiles=paste(c(ifile,list.files(paste(datadir,"/mcd09bias",sep=""),pattern=paste(sub("[.]tif","",basename(outfile3)),"_[0-9]*[.]tif",sep=""),full=T)),collapse=" ")
    system(paste("gdalbuildvrt -srcnodata 32767 -vrtnodata 32767 ",outfile," ",tfiles,sep="")) 
    system(paste("gdal_translate -a_ullr -180 90 180 -90 -a_nodata 32767 -co COMPRESS=LZW -co ZLEVEL=9 -co PREDICTOR=2 ",outfile," ",outfile2,sep=""))
    ## use pksetmask to set any values >100 (except the missing values) to 100
    ## these exist due to the reprojection to wgs84 from sinusoidal, there are a few pixels with values slightly
    ## greater than 100
    ## also convert to an unsigned 16-bit integer to allow adding a color table
    ict=ifelse(itype=="mean","data/colors16_mean.txt","data/colors16_sd.txt")
    system(paste("pksetmask -i ",outfile2," -m ",outfile2," -ot UInt16 ",
                 "--operator='>' --msknodata 20000 --nodata 65535  --operator='>' --msknodata 10000 --nodata 10000 --operator='<' --msknodata 0 --nodata 65535 ",
                 " -ct ",ict,"  -co COMPRESS=LZW -co PREDICTOR=2 -o ",outfile3))
    ## clean up temporary files
    file.remove(outfile,outfile2)
    writeLines(paste("#######################################             Finished ",outfile))
}


## check output
for(i in 1:nrow(df2)) {
    ifile=df2$path[i]
    outfile3=paste(datadir,"/mcd09ctif/",basename(ifile),sep="")
    print(ifile)
    system(paste("gdalinfo ",outfile3,"| grep 'Size is'"))
}
