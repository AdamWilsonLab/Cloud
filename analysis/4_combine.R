################################################################################
###  calculate monthly means of terra and aqua
source('analysis/setup.R')

beginCluster(12)

### assemble list of files to process
df=data.frame(path=list.files(paste(datadir,"/mcd09ctif",sep=""),full=T,pattern="M[Y|O].*[0-9]*[mean|sd].*tif$"),stringsAsFactors=F)
df[,c("sensor","month","type")]=do.call(rbind.data.frame,strsplit(basename(df$path),"_|[.]"))[,c(1,2,3)]

df2=unique(df[,c("month","type")])

overwrite=T

# Create combined (MOD+MYD) corrected mean CF
    foreach(i=1:nrow(df2), .options.multicore=list(preschedule=FALSE)) %dopar% {
        f=df$path[df$month==df2$month[i]&df$type==df2$type[i]]
        ## Define output and check if it already exists
        tmcd=paste(datadir,"/mcd09ctif/MCD09_",df2$type[i],"_",df2$month[i],"_uncompressed.tif",sep="")
        tmcd2=paste(datadir,"/mcd09ctif/MCD09_",df2$type[i],"_",df2$month[i],".tif",sep="")
        ## check if output already exists
        if(!overwrite&file.exists(tmcd2)){print(paste(tmcd2,"Exists, moving on..."));return(NULL)}
        if(overwrite&file.exists(tmcd2)){print(paste(tmcd2,"Exists, deleting it..."));file.remove(tmcd,tmcd2)}
        ## Take average between images
        ## switch NA values to 32768 to facilitate recasting to 8-bit below, otherwise they are confounded with 0 cloud values
        ops=paste(" -t_srs 'EPSG:4326' -multi -srcnodata 65535 -dstnodata 65535 -r bilinear -te -180 -90 180 90 -tr 0.008333333333333 -0.008333333333333",
            "-co BIGTIFF=YES  --config GDAL_CACHEMAX 20000 -wm 2000 -wo NUM_THREADS:10 -wo SOURCE_EXTRA=5")
        system(paste("gdalwarp -overwrite -r average -co COMPRESS=LZW -co ZLEVEL=9  ",ops," ",paste(f,collapse=" ")," ",tmcd))
        ## update metadata
        if(df2$type[i]=="mean")
            tags=c(paste("TIFFTAG_IMAGEDESCRIPTION='Monthly Cloud Frequency for 2000-2013 extracted from C5 MODIS MOD09GA and MYD09GA PGE11 internal cloud mask algorithm (embedded in state_1km bit 10).",
            "The daily cloud mask time series were summarized to mean cloud frequency (CF) by calculating the proportion of cloudy days.'"),
            "TIFFTAG_DOCUMENTNAME='Collection 5 MCD09 Cloud Frequency'",
            paste("TIFFTAG_DATETIME='2013",sprintf("%02d", i),"15'",sep=""),
              "TIFFTAG_ARTIST='Adam M. Wilson (adam.wilson@yale.edu)'")
        if(df2$type[i]=="sd")
        tags=c(paste("TIFFTAG_IMAGEDESCRIPTION='Standard Deviation of the Monthly Cloud Frequency for 2000-2013 extracted from C5 MODIS",
            " MOD09GA and MYD09GA PGE11 internal cloud mask algorithm (embedded in state_1km bit 10).",
            "The daily cloud mask time series were summarized to mean cloud frequency (CF) by calculating the proportion of cloudy days.'"),
            "TIFFTAG_DOCUMENTNAME='Collection 5 MCD09 SD of Cloud Frequency'",
            paste("TIFFTAG_DATETIME='2013",sprintf("%02d", i),"15'",sep=""),
              "TIFFTAG_ARTIST='Adam M. Wilson (adam.wilson@yale.edu)'")
        system(paste("/usr/local/src/gdal-1.10.0/swig/python/scripts/gdal_edit.py ",tmcd," ",paste("-mo ",tags,sep="",collapse=" "),sep=""))
        # create final fixed image
        system(paste("gdal_translate -co COMPRESS=LZW -co ZLEVEL=9 -co PREDICTOR=2 ",tmcd," ",tmcd2,sep=""))
        file.remove(tmcd)
        writeLines(paste("##########################################    Finished ",tmcd2))
    }


#################################################################################
###### Apply albedo mask and fill in missing data, convert to 8-bit compressed file, add colors and other details

f2=list.files(paste(datadir,"/mcd09ctif",sep=""),pattern=paste(".*MCD09_.*_[0-9].[.]tif$",sep=""),full=T)

## download albedo mask data from earth engine
download=F
hash="c5400b5cad"
if(download) system(paste("google docs get mask_",hash,"_2014066* ",datadir,"/mcd09ee_mask",sep=""))
mask=list.files(paste0(datadir,"/mcd09ee_mask"),pattern=paste0(".*",hash,".*tif$"),full=T)

## calculate what % of pixels were masked
system(paste0("gdalinfo -stats ",mask))
system(paste0("pkinfo  -src_min -128  -src_max 3   -hist -i ",mask," > data/out/mask_histogram.txt" ))
mhist=read.table("data/out/mask_histogram.txt")
mhist=mhist[mhist[,2]>0,]
mhist$p=100*mhist$V2/sum(mhist$V2)
100*sum(mhist$V2[mhist$V1>0])/sum(mhist$V2)

i=1


foreach(i=1:length(f2), .options.multicore=list(preschedule=FALSE)) %dopar% {
    file=f2[i]
    ## temporary files    
    tfile1=sub("[.]tif","_masked.tif",file)
    tfile1a=sub("[.]tif","_Wborder.tif",file)
    tfile1a2=sub("[.]tif","_Wborder2.tif",file)
    tfile1a3=sub("[.]tif","_Wborder3.tif",file)
    tfile1a4=sub("[.]tif","_Wborder4.tif",file)
    tfile1a5=sub("[.]tif","_Wborder5.tif",file)
    tfile1a6=sub("[.]tif","_Wborder6.tif",file)
    
    tfile1b=sub("[.]tif","_Eborder.tif",file)
    
    tfile2=sub("[.]tif","_filledrough.tif",file)
    tfile3=sub("[.]tif","_filledsmoothed.tif",file)
    tfile4=sub("[.]tif","_filled.tif",file)
    outfilevrt=sub("[.]tif",".vrt",file)
    outfile=paste("data/MCD09/",basename(file),sep="")
    ## mask the image
    system(paste0("pksetmask -i ",file," -m ",mask," --operator='>' --msknodata 0 --nodata 32767 -ot UInt16 -o ",tfile1))
    ## Subset edges out for special treatment to avoid errors in gdal_fillnodata
    system(paste0("gdal_translate -a_nodata 32767 -ot UInt16 -srcwin 43198 0 2 21600 ",file," ",tfile1a))
    system(paste0("gdal_translate -a_nodata 32767 -ot UInt16 -srcwin 43195 0 5 21600 ",file," ",tfile1a2))
    
    #system(paste0("gdal_translate -srcwin 0 0 2 21600 ",file," ",tfile1b))
    ## replace any missing values (65535) with 32767 for filling
    system(paste0("pkfilter -f max -dx 3 -dy 1 -i ",tfile1a,"  -o ",tfile1a3)) #",tfile1,"
    system(paste0("pkfilter -f mean -dx 5 -dy 1 -i ",tfile1a2,"   -o ",tfile1a4))
    system(paste0("gdal_translate -srcwin 3 0 2 21600 ",tfile1a4," ",tfile1a5))
    system(paste0("pksetmask -i ",tfile1a5," -m ",tfile1a3," --operator='>' --msknodata 10000 --nodata 32767 -o ",tfile1a6))
    system(paste0("gdal_edit.py -a_nodata 32767 ",tfile1a6))
    
    ## merge this back into the masked file
    system(paste0("gdalwarp ",tfile1a6," ",tfile1))    
    ## fill no data
    system(paste0("gdal_fillnodata.py -md 1000 ",tfile1," ",tfile2))
    #system(paste0("pkfillnodata -co 'BIGTIFF=TRUE' -co 'COMPRESS=LZW' -d 1000 -i ",file," -m ",tfile1," -o ",tfile2))
    
    ## smooth the filled data
    system(paste0("pkfilter -f smooth -dx 11 -dy 11 -nodata 32767  -i ",tfile2," -o ",tfile3))
    ## used the smoothed data rather than the original, would prefer to use the -si option on gdal_fillnodata, but it's broken
    system(paste0("gdal_merge.py -o ",tfile4," -co COMPRESS=LZW -co PREDICTOR=2 -n 32767 -init 32767 ",tfile3," ",tfile1))    
#-scale 0 10000 0 100 
    system(paste("gdal_translate  -of VRT ",tfile4," ",outfilevrt)) 
    ## add color table for 8-bit data
    vrt=scan(outfilevrt,what="char")
    hd=c("<ColorInterp>Palette</ColorInterp>","<ColorTable>")
    ft="</ColorTable>"
    colR=colorRampPalette(c("#08306b","#0d57a1","#2878b8","#4997c9","#72b2d7","#a2cbe2","#c7dcef","#deebf7","#f7fbff"))
    cols=data.frame(t(col2rgb(colR(10005))))
    ct=paste("<Entry c1=\"",cols$red,"\" c2=\"",cols$green,"\" c3=\"",cols$blue,"\" c4=\"32767\"/>")
    cti=grep("ColorInterp",vrt)  # get index of current color table
    vrt2=c(vrt[1:(cti-1)],hd,ct,ft,vrt[(cti+1):length(vrt)])
    ## update missing data flag following http://lists.osgeo.org/pipermail/gdal-dev/2010-February/023541.html
    #csi=grep("<ComplexSource>",vrt2)  # get index of current color table
    #vrt2=c(vrt2[1:csi],"<NODATA>32767</NODATA>",vrt2[(csi+1):length(vrt2)])
    write.table(vrt2,file=outfilevrt,col.names=F,row.names=F,quote=F)              
    tags=c(paste("TIFFTAG_IMAGEDESCRIPTION='Monthly Cloud Frequency for 2000-2013 extracted from C5 MODIS M*D09GA PGE11 internal cloud mask algorithm (embedded in state_1km bit 10).",
                 "The daily cloud mask time series were summarized to mean cloud frequency (CF) by calculating the proportion of cloudy days.'"),
           "TIFFTAG_DOCUMENTNAME='Collection 5 Cloud Frequency'",
           paste("TIFFTAG_DATETIME='2014'",sep=""),
           "TIFFTAG_ARTIST='Adam M. Wilson (adam.wilson@yale.edu)'")
    system(paste("gdal_translate -a_nodata 65535 -ot UInt16  -co COMPRESS=LZW -co PREDICTOR=2 ",paste("-mo ",tags,sep="",collapse=" ")," ",outfilevrt," ",outfile))
  if(file.exists(outfile)) print(paste("Finished ",outfile))
    ## clean up  - keep the 16 bit filled version for additional calculations
    file.remove(tfile1,tfile1a,tfile1a2,tfile1a3,tfile1a4,tfile1a5,tfile2,tfile3)
}

########################################################################
f3=list.files(paste("data/MCD09/",sep=""),pattern=paste("MCD09_mean_[0-9].[.]tif$",sep=""),full=T)
f3sd=list.files(paste("data/MCD09/",sep=""),pattern=paste(".*MCD09_sd_[0-9].[.]tif$",sep=""),full=T)

##############################################################################
### create a temporary copy without a color table for uploading to earth engine
dir.create("data/MCD09_EarthEngineUpload")

foreach(f=c(f3,f3sd))%dopar% {
  tvrt=paste0("data/MCD09_EarthEngineUpload/",sub("tif","vrt",basename(f)))
  tout=paste0("data/MCD09_EarthEngineUpload/",basename(f))
  system(paste0("gdal_translate -ot UInt16 -of vrt ",f," ",tvrt))
  yvrt=scan(tvrt,what="char",)
  cti=grep("ColorTable",yvrt)
  yvrt2=c(yvrt[1:(cti[1]-1)],"<ColorInterp>Grey</ColorInterp>",yvrt[(cti[2]+1):length(yvrt)])
  write.table(yvrt2,file=tvrt,col.names=F,row.names=F,quote=F)              
  system(paste("gdal_translate -a_nodata 65535 -ot UInt16  -co COMPRESS=LZW -co PREDICTOR=2 ",tvrt," ",tout))  
  file.remove(tvrt)
}
  

dmean=stack(as.list(f3))
NAvalue(dmean)=65535
#dmean=setZ(dmean,z=as.Date(paste0("2014-",1:12,"-15"))-as.Date("2000-01-01"))

dsd=stack(as.list(f3sd))
NAvalue(dsd)=65535

################
### calculate inter vs. intra annual variability

## Function to calculate standard deviation and round it to nearest integer
Rsd=function(x) calc(x,function(x) {
  if(all(is.na(x))) return(NA)
  sd(x,na.rm=T)
})


Rmean=function(x) calc(x,function(x) {
  mean(x,na.rm=T)
})


dintra=clusterR(dmean,Rsd,m=4,file="data/MCD09_deriv/intra_uncompressed.tif",
                #options=c("COMPRESS=LZW","PREDICTOR=2"),
                overwrite=T,datatype='INT2U',NAflag=65535)
system(paste("gdal_translate -ot UInt16 -co COMPRESS=DEFLATE -co ZLEVEL=9 -co BIGTIFF=YES -co COMPRESS=LZW -co PREDICTOR=2",
" data/MCD09_deriv/intra_uncompressed.tif data/MCD09_deriv/intra.tif"))
file.remove("data/MCD09_deriv/intra_uncompressed.tif")


dinter=clusterR(dsd,Rmean,m=4,file="data/MCD09_deriv/inter_uncompressed.tif",
                #options=c("COMPRESS=LZW","PREDICTOR=2"),
                overwrite=T,datatype='INT2U',NAflag=65535)
system(paste("gdal_translate -ot UInt16 -co COMPRESS=DEFLATE -co ZLEVEL=9 -co BIGTIFF=YES -co COMPRESS=LZW -co PREDICTOR=2",
       " data/MCD09_deriv/inter_uncompressed.tif data/MCD09_deriv/inter.tif"))
file.remove("data/MCD09_deriv/inter_uncompressed.tif")


## Overall annual mean
dmeanannual=clusterR(dmean,Rmean,m=4,file="data/MCD09_deriv/meanannual_uncompressed.tif",
#                     options=c("COMPRESS=LZW","PREDICTOR=2"),
                     overwrite=T,datatype='INT2U',NAflag=65535)
system(paste("gdal_translate -ot UInt16 -co COMPRESS=DEFLATE -co ZLEVEL=9 -co BIGTIFF=YES -co COMPRESS=LZW -co PREDICTOR=2",
             " data/MCD09_deriv/meanannual_uncompressed.tif data/MCD09_deriv/meanannual.tif"))
file.remove("data/MCD09_deriv/meanannual_uncompressed.tif")


#file.copy("data/MCD09_deriv/meanannual.tif","data/MCD09_EarthEngineUpload/meanannual.tif")
#file.copy("data/MCD09_deriv/inter.tif","data/MCD09_EarthEngineUpload/inter.tif")
#file.copy("data/MCD09_deriv/intra.tif","data/MCD09_EarthEngineUpload/intra.tif")

#################################################
### Calculate Markham's Seasonality

#fseasconc(dmean)
#stheta=calc(dmean,seastheta)

fseasconc=function(x) calc(x,seasconc,na.rm=F)
fseastheta=function(x) calc(x,seastheta,na.rm=F)

## calculate seasonality
sconc=clusterR(dmean,fun=fseasconc,export="seasconc",overwrite=T,filename="data/MCD09_deriv/seasconc.tif",NAflag=65535,datatype="INT2U")
stheta=clusterR(dmean,fun=fseastheta,export="seastheta",overwrite=T,filename="data/MCD09_deriv/seastheta.tif",NAflag=65535,datatype="INT2U")


### generate color key
seas=stack("data/MCD09_deriv/seasconc.tif","data/MCD09_deriv/seastheta.tif")
gain(seas)=.1
names(seas)=c("conc","theta")
=======
sconc=clusterR(dmean,
               fun=fseasconc,
               export="seasconc",
               overwrite=T,
               filename="data/MCD09_deriv/seasconc.tif",
               NAflag=65535,
               datatype="INT2U")

stheta=clusterR(dmean,
                fun=fseastheta,
                export="seastheta",
                overwrite=T,
                filename="data/MCD09_deriv/seastheta.tif",
                NAflag=65535,
                datatype="INT2U")


### generate color key
seas=stack("data/MCD09_deriv/seasconc.tif",
           "data/MCD09_deriv/seastheta.tif")
gain(seas)=.1
names(seas)=c("conc",
              "theta")
NAvalue(seas)=65535

## extract unique values from seas to develop color table
## takes ~30 minutes
parUnique=function(x){
  ## break up raster into chunks for more efficient unique value search
  blks=blockSize(x, n=nlayers(seas), minblocks=200, minrows=1)
  uv=foreach(i=1:blks$n,.combine=rbind,.inorder=F)%dopar%{
    # extract values for this chunk
    tvals=round(getValues(x, blks$row[i], blks$nrows[i])/10)
    ## convert to single vector for faster unique()
    tvals2=paste(tvals[,1],tvals[,2],sep="_")
    ## split back to separate values
    tuvals=do.call(rbind,strsplit(unique(tvals2),split="_"))
    rcols=data.frame(conc=as.numeric(tuvals[,1]),theta=as.numeric(tuvals[,2]))  
    return(rcols)
  }
  ## find overall unique values
  uvals=unique(uv)
}

fuvals="data/out/uniqueseasonality.csv"
if(!file.exists(fuvals)){
  uvals=parUnique(seas)
  write.csv(uvals,file=fuvals,row.names=F)
}

## read in unique colors
uvals=read.csv(fuvals)


## Generate color table
summary(uvals)
concs=seq(0,76,len=150)
thetas=seq(0,360,len=360)
col=expand.grid(conc=concs,theta=thetas)
col=cbind.data.frame(col,t(col2rgb(as.character(cut(col$theta,breaks=1000,labels=rainbow(1000))))))
col=cbind.data.frame(col,t(rgb2hsv(r=col$red,g=col$green,b=col$blue,maxColorValue=255)))
col$id=as.integer(1:nrow(col))
prot=-180
col$x=col$conc*cos((pi/180)*(-col$theta+prot))
col$y=col$conc*sin((pi/180)*(-col$theta+prot))
## adjust saturation and value to bring low concentration values down
cbreak=10
col$s=as.numeric(as.character(cut(col$conc,breaks=c(0,seq(1,cbreak*1.5,len=51),seq(cbreak*1.5+1,max(concs),len=50)),labels=c(0,seq(0.3,1,len=100)))))
col$v=as.numeric(as.character(cut(col$conc,breaks=c(0,seq(1,cbreak*1.5,len=51),seq(cbreak*1.5+1,max(concs),len=50)),labels=c(0,seq(0.01,1,len=100)))))
col$s[is.na(col$s)]=0;col$v[is.na(col$v)]=0
col$val=hsv(h=col$h,s=col$s,v=col$v)
## rewrite RGB values to update saturation and value
col[c("r","g","b","alpha")]=t(col2rgb(col$val,alpha=T))
col$rgbcol=rgb(red=col$r,blue=col$b,green=col$g,max=255)

col$exists=paste(round(col$conc),round(col$theta))%in%paste(uvals$conc,uvals$theta)

## Write the color key table
write.csv(col,row.names=F,file="data/MCD09_deriv/coltable.csv")

## create new raster referencing this color table
tcol=as(dcast(col,conc~theta,value.var="id",df=F,fill=NA),"matrix")
dimnames(tcol)=NULL
## function to look up color table value for each pixel
iseas=function(x,...){
  calc(x,function(v,...) {
    if(any(is.na(v))) return(NA)
    tcol[which.min(abs(concs-v[1])),which.min(abs(thetas-v[2]))]
    })
}
## run the function to find the color ID for each pixel
## must not use compression flag (COMPRESS=LZW) or geotif export will fail
seas2 <- clusterR(seas, iseas, export=c('concs','thetas','tcol'),
                  dataType="INT2U",filename="data/MCD09_deriv/seas_ind.tif",
                  NAflag=65534,na.rm=T,overwrite=T)

endCluster()

## update color table via VRT
seasvrt="data/MCD09_deriv/seas_vis.vrt"
## write a VRT
system(paste("gdal_translate -of VRT -a_nodata 65534 data/MCD09_deriv/seas_ind.tif ",seasvrt)) 
## read in the VRT and update the color information using the color table
vrt=scan(seasvrt,what="char")
hd=c("<ColorInterp>Palette</ColorInterp>","<ColorTable>")
ft="</ColorTable>"
ct=paste("<Entry c1=\"",col$r,"\" c2=\"",col$g,"\" c3=\"",col$b,"\" c4=\"255\"/>")
cti=grep("ColorInterp",vrt)  # get index of current color table
vrt2=c(vrt[1:(cti-1)],hd,ct,ft,vrt[(cti+1):length(vrt)])
write.table(vrt2,file=seasvrt,col.names=F,row.names=F,quote=F)              
## convert back to geotif and compress
system(paste("gdal_translate -co COMPRESS=LZW -a_nodata 65535 -co PREDICTOR=2 -of GTIFF -ot UInt16 ",seasvrt,"  data/MCD09_deriv/seas_visct.tif")) 

## expand to RGB for EarthEngine
system(paste("gdal_translate -a_nodata 255 -expand rgba -co COMPRESS=LZW -co PREDICTOR=2 -of GTIFF -ot Byte ",seasvrt,"  data/MCD09_deriv/seas_rgb.tif")) 

#file.copy("data/MCD09_deriv/seas_rgb.tif","data/MCD09_EarthEngineUpload/seas_rgb.tif")
