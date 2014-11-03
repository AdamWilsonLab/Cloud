
# Load the libraries and set working directory
source("analysis/setup.R")



####################################################################
###  Map of hotspots
  registerDoMC(12)
  qs=c(0,0.01,0.025,0.05,0.1,0.25,0.75,0.9,0.95,.975,0.99,1)
  ## list of products to process
  qprods=list("data/MCD09_deriv/meanannual.tif",
              "data/MCD09_deriv/inter.tif",
              "data/MCD09_deriv/intra.tif",
              "data/MCD09_deriv/mean_1deg_sd.tif")
  
  ## get elevation data for land mask
  land=raster("/mnt/data/jetzlab/Data/environ/global/worldclim/alt.bil")
  NAvalue(land)=-9999
  land=crop(land,extent(greg$xlim,greg$ylim))

  system("gdal_translate -co COMPRESS=LZW /mnt/data/jetzlab/Data/environ/global/worldclim/alt.bil data/worldclim_alt.tif")
  
  ## make table of quantiles
  if(!file.exists("output/derivedquantiles.csv")){
    qsr=foreach(i=1:length(qprods),.combine=rbind.data.frame)%dopar%{  
      tfile=paste0(tempdir(),"/",sub(".tif",".vrt",basename(qprods[[i]])))     
      ## crop to region of interest (removing antarctica)
      system(paste("gdalwarp -of vrt -te ",
                   paste(c(greg$xlim[1],greg$ylim[1],greg$xlim[2],greg$ylim[2]),collapse=" "),
                   qprods[[i]],tfile))
      ## crop to land
      trast=raster(tfile)
      trast=mask(trast,land)
      NAvalue(trast)=0
      tq=quantile(trast, probs=qs,ncells=NULL,na.rm=T,names=T)
      rmr(trast)
      file.remove(tfile)
      ## calculate quantiles
      cbind.data.frame(name=basename(qprods[[i]]),t(tq))
    }
    colnames(qsr)=c("name",paste0("Q",qs))
    
  write.csv(qsr,"output/derivedquantiles.csv",row.names=F)
  }

  ## read in quantile tables
  qsr=read.csv("output/derivedquantiles.csv")
  
  qsrc=data.frame(id=  c(0, 1,  2,   3,       4:255),
                  R=   c(0, 0,  192, 255, rep(0,252)),
                  G=   c(0, 0,  192, 0  , rep(0,252)),
                  B=   c(0, 255,192, 0  , rep(0,252)),
                  ALFA=c(0, 255,255, 255, rep(0,252)))
  write.table(qsrc,"output/derivedquantiles_colors.csv",quote=F,col.names=F,row.names=F)
  
  
  
  ## reclass each raster to low-medium-high
    foreach(i=1:length(qprods))%dopar%{  
#      tfile=paste0(tempdir(),"/",sub(".tif",".vrt",basename(qprods[[i]])))
      outfile=paste0("data/MCD09_deriv/",sub(".tif","_hotspots.tif",basename(qprods[[i]])))     
      
      if(file.exists(outfile)) return(NA)
      tq=qsr[qsr$name==basename(qprods[[i]]),]#sub(".vrt",".tif",basename(tfile)),]        
      ## crop to region of interest (removing antarctica)
#      system(paste("gdalwarp -of vrt -te ",
#                   paste(c(greg$xlim[1],greg$ylim[1],greg$xlim[2],greg$ylim[2]),collapse=" "),
#                   qprods[[i]],tfile))
      system(paste0("pksetmask -ot Byte -i ",qprods[[i]],
                    " -m data/worldclim_alt.tif --operator='<' --msknodata -1000 --nodata 0 ",                    
                    " -m ",qprods[[i]]," --operator='>' --msknodata ",tq$Q1," --nodata 0 ", # missing data
                    " -m ",qprods[[i]]," --operator='<' --msknodata ",tq$Q0.025," --nodata 1 ", # low 
                    " -m ",qprods[[i]]," --operator='<' --msknodata ",tq$Q0.975," --nodata 2 ", # middle 
                    " -m ",qprods[[i]]," --operator='>' --msknodata ",tq$Q0.975," --nodata 3 ", # top 
                    " -ct output/derivedquantiles_colors.csv ",
                    " -o ",outfile))
    }
    
  }
  
  
 
## find which combinations actually exist
  hs=stack(list.files(paste0("data/MCD09_deriv/"),pattern="hotspots",full=T))

## build category table
cnames=factor(c(1,2,3),ordered=T) #,labels=c("low","mid","high")
cgrid=expand.grid(inter=cnames,intra=cnames,spatial=cnames,mean=cnames)

# function to assign each gric cell to row in cgrid
getid=function(x) which(apply(cgrid, 1, function(y) all(y == x)))

calc(hs,fun=getid,file="data/MCD09_deriv/combined_hotspots.tif")

  
}