################################################################################
###  calculate monthly means of terra and aqua
source('analysis/setup.R')

library(raster)
library(foreach)
library(multicore)
library(doMC)
registerDoMC(20)

beginCluster(20)

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
###### convert to 8-bit compressed file, add colors and other details


f2=list.files(paste(datadir,"/mcd09ctif",sep=""),pattern=paste(".*MCD09_.*_[0-9].[.]tif$",sep=""),full=T)


foreach(i=1:length(f2), .options.multicore=list(preschedule=FALSE)) %dopar% {
    file=f2[i]
    outfilevrt=sub("[.]tif",".vrt",file)
    outfile=paste("data/mcd09tif/",basename(file),sep="")
    ## rescale to 0-100 using a VRT
    system(paste("gdal_translate  -scale 0 10000 0 100 -of VRT ",file," ",outfilevrt)) 
    ## add color table for 8-bit data
    vrt=scan(outfilevrt,what="char")
    hd=c("<ColorInterp>Palette</ColorInterp>","<ColorTable>")
    ft="</ColorTable>"
    colR=colorRampPalette(c("#08306b","#0d57a1","#2878b8","#4997c9","#72b2d7","#a2cbe2","#c7dcef","#deebf7","#f7fbff"))
    cols=data.frame(t(col2rgb(colR(105))))
    ct=paste("<Entry c1=\"",cols$red,"\" c2=\"",cols$green,"\" c3=\"",cols$blue,"\" c4=\"255\"/>")
    cti=grep("ColorInterp",vrt)  # get index of current color table
    vrt2=c(vrt[1:(cti-1)],hd,ct,ft,vrt[(cti+1):length(vrt)])
    ## update missing data flag following http://lists.osgeo.org/pipermail/gdal-dev/2010-February/023541.html
    csi=grep("<ComplexSource>",vrt2)  # get index of current color table
    vrt2=c(vrt2[1:csi],"<NODATA>327</NODATA>",vrt2[(csi+1):length(vrt2)])
    write.table(vrt2,file=outfilevrt,col.names=F,row.names=F,quote=F)              
    tags=c(paste("TIFFTAG_IMAGEDESCRIPTION='Monthly Cloud Frequency for 2000-2013 extracted from C5 MODIS M*D09GA PGE11 internal cloud mask algorithm (embedded in state_1km bit 10).",
        "The daily cloud mask time series were summarized to mean cloud frequency (CF) by calculating the proportion of cloudy days. ",
        "Band Descriptions: 1) Mean Monthly Cloud Frequency'"),
        "TIFFTAG_DOCUMENTNAME='Collection 5 Cloud Frequency'",
        paste("TIFFTAG_DATETIME='2014'",sep=""),
        "TIFFTAG_ARTIST='Adam M. Wilson (adam.wilson@yale.edu)'")
    system(paste("gdal_translate -a_nodata 255 -ot Byte  -co COMPRESS=LZW -co PREDICTOR=2 ",paste("-mo ",tags,sep="",collapse=" ")," ",outfilevrt," ",outfile))

}


################
### calculate inter vs. intra annual variability
f3=list.files(paste(datadir,"/mcd09ctif/",sep=""),pattern=paste(".*MCD09_mean_[0-9].[.]tif$",sep=""),full=T)
f3sd=list.files(paste(datadir,"/mcd09ctif/",sep=""),pattern=paste(".*MCD09_sd_[0-9].[.]tif$",sep=""),full=T)

dmean=stack(as.list(f3))
dsd=stack(as.list(f3sd))


## Function to calculate standard deviation and round it to nearest integer
Rsd=function(x) calc(x,function(x) {
  sd=sd(x,na.rm=T)
  if(is.na(sd)) sd=0
  return(round(sd/100))
})

Rmean=function(x) calc(x,function(x) {
  return(round(mean(x,na.rm=T)/100))
})

dintra=clusterR(dmean,Rsd,na.rm=T,file="data/MCD09_deriv/intra.tif",options=c("COMPRESS=LZW","PREDICTOR=2"),overwrite=T,dataType='INT1U',NAflag=255)
dinter=clusterR(dsd,Rmean,na.rm=T,file="data/MCD09_deriv/inter.tif",options=c("COMPRESS=LZW","PREDICTOR=2"),overwrite=T,dataType='INT1U',NAflag=255)

## Overall annual mean
dmeanannual=clusterR(dmean,Rmean,na.rm=T,file="data/MCD09_deriv/meanannual.tif",options=c("COMPRESS=LZW","PREDICTOR=2"),overwrite=T,dataType='INT1U',NAflag=255)


#################################################
### Calculate Markham's Seasonality
#tdmean=crop(dmean,extent(c(17,30,-35,-28)))


#seasconc(dmean[400000],return.Pc=T,return.thetat=T)

## calculate seasonality
seas=calc(dmean,seasconc ,overwrite=T,
                    options=c("COMPRESS=LZW","PREDICTOR=2"),
                    filename="data/MCD09_deriv/seas_conc.tif",NAflag=65535,datatype="INT2U")



### generate color key
seas=stack("data/MCD09_deriv/seas_conc.tif")
gain(seas)=.1
names(seas)=c("conc","theta")

#r="Venezuela"
#seas=crop(seas,regs[[r]])


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

## Set polar rotation for polar plot of color values
prot=-180

## Generate color table
summary(uvals)
concs=seq(0,76,len=150)
thetas=seq(0,360,len=360)
col=expand.grid(conc=concs,theta=thetas)
col=cbind.data.frame(col,t(col2rgb(as.character(cut(col$theta,breaks=1000,labels=rainbow(1000))))))
col=cbind.data.frame(col,t(rgb2hsv(r=col$red,g=col$green,b=col$blue,maxColorValue=255)))
col$id=as.integer(1:nrow(col))
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
tcol=as(cast(col,conc~theta,value="id",df=F,fill=NA),"matrix")
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
system(paste("gdal_translate -co COMPRESS=LZW -a_nodata 65534 -co PREDICTOR=2 -of GTIFF -ot UInt16 ",seasvrt,"  data/MCD09_deriv/seas_visct.tif")) 

## expand to RGB for EarthEngine
system(paste("gdal_translate -a_nodata 255 -expand rgba -co COMPRESS=LZW -co PREDICTOR=2 -of GTIFF -ot Byte ",seasvrt,"  data/MCD09_deriv/seas_rgb.tif")) 


pCirc=function(r,n=100) {
  angs=seq(0,360,len=100)
  #lapply(r,function(r) {
      cbind(
      y=r*cos((pi/180)*(-angs)),
      x=r*sin((pi/180)*(-angs))
      )
  #})
}

mangle=(-seq(30,360,30)+prot+15)*(pi/180)
rmons=c(1:3,9:12)
lmons=4:8
mon=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec") # month names
ladj=65
lims=c(-75,75)

png(width=1000,height=1000,pointsize=32,file="manuscript/figures/SeasKey_%0d.png")
### Create the color key
xyplot(conc~theta,col=col$val[col$exists],data=col[col$exists,],pch=16,cex=1,
       scales=list(
         x=list(labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
                at=seq(15,360,30))),
       xlab="Month",ylab="Concentration (%)")

xyplot(y~x,col=col$val[col$exists],data=col[col$exists,],pch=16,cex=1.5,
       xlab="",ylab="",scales=list(draw=F),ylim=lims,xlim=lims,asp=1,
       par.settings = list(axis.line = list(col = "transparent")))+
  layer(panel.polygon(pCirc(r=20,n=100),border="grey"))+
  layer(panel.polygon(pCirc(r=40,n=100),border="grey"))+
  layer(panel.polygon(pCirc(r=60,n=100),border="grey"))+
  layer(panel.segments(0,0,60*cos(mangle-(15*pi/180)),60*sin(mangle-(15*pi/180)),col="grey"))+ #draw angles
  layer(panel.text(x=ladj*cos(mangle[lmons]),y=ladj*sin(mangle[lmons]),mon[lmons],pos=4,cex=3,offset=0,srt=(mangle[lmons])*180/pi))+ #add left months))
  layer(panel.text(x=ladj*cos(mangle[rmons]),y=ladj*sin(mangle[rmons]),mon[rmons],pos=2,cex=3,offset=0,srt=((mangle[rmons])*180/pi)-180))+ # add right months
  layer(panel.text(x=-2,y=c(-20,-40,-60),c(20,40,60),pos=4,cex=3)) #add scale text

dev.off()


## standard graphics plotting
#plot(col$x,col$y,xaxt="n",yaxt="n",xlab="",ylab="",col=col$val,pch=16,cex=1.5,new=F,asp=1,bty="n",xlim=c(-16,16),ylim=c(-40,40))
#plot(col$x[col$exists],col$y[col$exists],xaxt="n",yaxt="n",xlab="",ylab="",col=col$val[col$exists],pch=16,cex=2,new=F,asp=1,bty="n",xlim=c(-16,16),ylim=c(-45,45))
#n=10; circ=seq(0,30,n) #define the circle diameters
#symbols(rep(0,length(circ)),rep(0,length(circ)),circles=circ,add=T,fg="grey",lwd=2,inches=F) #draw circles
#segments(0,0,30*cos(mangle),30*sin(mangle),col="grey") #draw angles
#text(x=-2,y=-circ[2:n],circ[2:n],pos=4,cex=1.5) #add scale text
#mon=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec") # month names
#ladj=35
#for(i in c(1:3,9:12)) text(x=ladj*cos(-mangle[i]),y=ladj*sin(-mangle[i]),mon[i],pos=4,cex=1.8,offset=0,srt=-mangle[i]*180/pi) #add left months
#for(i in 4:8) text(x=ladj*cos(-mangle[i]),y=ladj*sin(-mangle[i]),mon[i],pos=2,cex=1.8,offset=0,srt=(-mangle[i]*180/pi)-180) # add right months

########################################################################################
#### stuff below here is old junk.....       
        
        ##  convert to netcdf, subset to mean/sd bands
        trans_ops=paste(" -co COMPRESS=DEFLATE -a_nodata 255 -stats -co FORMAT=NC4 -co ZLEVEL=9 -b 4 -b 2")
        system(paste("gdal_translate -of netCDF ",trans_ops," ",ttif2," ",ncfile))
        ## file.remove(temptffile)
        system(paste("ncecat -O -u time ",ncfile," ",ncfile,sep=""))
        ## create temporary nc file with time information to append to CF data
        cat(paste("
    netcdf time {
      dimensions:
        time = 1 ;
      variables:
        int time(time) ;
      time:units = \"days since 2000-01-01 ",ifelse(s=="MOD09","10:30:00","13:30:00"),"\" ;
      time:calendar = \"gregorian\";
      time:long_name = \"time of observation\"; 
    data:
      time=",as.integer(date-as.Date("2000-01-01")),";
    }"),file=paste(tempdir(),"/",date,"_time.cdl",sep=""))
        system(paste("ncgen -o ",tempdir(),"/",date,"_time.nc ",tempdir(),"/",date,"_time.cdl",sep=""))
        ## add time dimension to ncfile and compress (deflate)
        system(paste("ncks --fl_fmt=netcdf4 -L 9 -A ",tempdir(),"/",date,"_time.nc ",ncfile,sep=""))
        ## add other attributes
        system(paste("ncrename -v Band1,CF ",ncfile,sep=""))
        system(paste("ncrename -v Band2,CFsd ",ncfile,sep=""))
        ## build correction factor explanation
        system(paste("ncatted ",
                     ## CF Mean
                     " -a units,CF,o,c,\"%\" ",
                     " -a valid_range,CF,o,ub,\"0,100\" ",
                                        #               " -a scale_factor,CF,o,b,\"0.1\" ",
                     " -a _FillValue,CF,o,ub,\"255\" ",
                     " -a missing_value,CF,o,ub,\"255\" ",
                     " -a long_name,CF,o,c,\"Cloud Frequency (%)\" ",
                     " -a correction_factor_description,CF,o,c,\"To account for variable observation frequency, CF in each pixel was adjusted by the proportion of days with at least one MODIS observation\" ",
                     " -a correction_factor,CF,o,f,\"",round(modbeta1,4),"\" ",
                     " -a NETCDF_VARNAME,CF,o,c,\"Cloud Frequency (%)\" ",
                     ## CF Standard Deviation
                     " -a units,CFsd,o,c,\"SD\" ",
                     " -a valid_range,CFsd,o,ub,\"0,200\" ",
                     " -a scale_factor,CFsd,o,f,\"0.25\" ",
                     " -a _FillValue,CFsd,o,ub,\"255\" ",
                     " -a missing_value,CFsd,o,ub,\"255\" ",
                     " -a long_name,CFsd,o,c,\"Cloud Frequency (%) Intra-month (2000-2013) Standard Deviation\" ",
                     " -a NETCDF_VARNAME,CFsd,o,c,\"Cloud Frequency (%) Intra-month (2000-2013) Standard Deviation\" ",
                     ## global
                     " -a title,global,o,c,\"Cloud Climatology from MOD09 Cloud Mask\" ",
                     " -a institution,global,o,c,\"Jetz Lab, EEB, Yale University, New Haven, CT\" ",
                     " -a source,global,o,c,\"Derived from MOD09GA Daily Data\" ",
                     " -a comment,global,o,c,\"Developed by Adam M. Wilson (adam.wilson@yale.edu / http://adamwilson.us)\" ",
                     ncfile,sep=""))

        ## add the fillvalue attribute back (without changing the actual values)
                                        #system(paste("ncatted -a _FillValue,CF,o,b,-32768 ",ncfile,sep=""))
        
        if(as.numeric(system(paste("cdo -s ntime ",ncfile),intern=T))<1) {
            print(paste(ncfile," has no time, deleting"))
            file.remove(ncfile)
        }
        print(paste(basename(ncfile)," Finished"))
        
        
    }

if(ramdisk) {
    ## unmount the ram disk
    system(paste("sudo umount ",tmpfs)
}


### merge all the tiles to a single global composite
#system(paste("ncdump -h ",list.files(tempdir,pattern="mod09.*.nc$",full=T)[10]))
file.remove("tmp/mod09_2000-01-15.nc")
system(paste("cdo -O  mergetime -setrtomiss,-32768,-1 ",paste(list.files(tempdir,pattern="mod09.*.nc$",full=T),collapse=" ")," data/cloud_monthly.nc"))

#  Overall mean
system(paste("cdo -O  timmean data/cloud_monthly.nc  data/cloud_mean.nc"))

## Seasonal Means
system(paste("cdo  -f nc4c -O -yseasmean data/cloud_monthly.nc data/cloud_yseasmean.nc"))
system(paste("cdo  -f nc4c -O -yseasstd data/cloud_monthly.nc data/cloud_yseasstd.nc"))


## standard deviation of mean monthly values give intra-annual variability
system("cdo -f nc4c -z zip -chname,CF,CFsd -timstd data/cloud_ymonmean.nc data/cloud_std_intra.nc")
## mean of monthly standard deviations give inter-annual variability 
system("cdo -f nc4c -z zip -chname,CF,CFsd -timmean data/cloud_ymonstd.nc data/cloud_std_inter.nc")




## Daily animations
regs=list(
    Venezuela=extent(c(-69,-59,0,7)),
    Cascades=extent(c(-122.8,-118,44.9,47)),
    Hawaii=extent(c(-156.5,-154,18.75,20.5)),
    Boliva=extent(c(-71,-63,-20,-15)),
    CFR=extent(c(17.75,22.5,-34.8,-32.6)),
    Madagascar=extent(c(46,52,-17,-12))
    )

r=1

system(paste("cdo  -f nc4c -O inttime,2012-01-15,12:00:00,7day  -sellonlatbox,",
             paste(regs[[r]]@xmin,regs[[r]]@xmax,regs[[r]]@ymin,regs[[r]]@ymax,sep=","),
             "  data/cloud_monthly.nc data/daily_",names(regs[r]),".nc",sep=""))




