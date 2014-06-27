
## libraries
library(rasterVis)
library(latticeExtra)
library(xtable)
library(texreg)
library(reshape2)
library(caTools)
library(rgeos)
library(raster)
library(plyr)
library(knitr)
require(knitcitations)
## read in global coasts for nice plotting
library(maptools)
library(rgdal)
library(coda)
#library(spgrass6)

## install rasterAutocorr
#library(devtools) 
#install_github("adammwilson/rasterAutocorr")
library(rasterAutocorr)

## register parallel backend
library(doMC)
registerDoMC(18)

rasterOptions(progress="text",maxmemory=1e12)

## increase the amount of memory for gdal to speed up processing
Sys.setenv(GDAL_CACHEMAX=5000,CPL_LOG_ERRORS="ON")

## set working directory for different machines
if( Sys.info()["nodename"]=="repens"){
  setwd("/Users/adamw/repos/Cloud")
  # temporary files will be written here:
  #  datadir="/mnt/data2/projects/cloud/"
}

if( Sys.info()["nodename"]=="litoria"){  
  setwd("/media/data/Cloud")
  # temporary files will be written here:
  datadir="/mnt/data2/projects/cloud/"
}

## color ramps
colR=colorRampPalette(c("#08306b","#0d57a1","#2878b8","#4997c9","#72b2d7","#a2cbe2","#c7dcef","#deebf7","#f7fbff"))
bgr=colorRampPalette(c("#0000ff","#00ff00","#ff0000"))
bgyrp=colorRampPalette(c("blue","darkgreen","goldenrod","red","purple"))

bgr=function(x,n=100,br=0,c1=c("darkblue","blue","grey"),c2=c("grey","red","purple")){
  at=unique(c(seq(min(x,na.rm=T),max(x,na.rm=T),len=n)))
  bg=colorRampPalette(c1)
  gr=colorRampPalette(c2)
  return(list(at=at,col=c(bg(sum(at<br)),gr(sum(at>=br)))))
}


## Set polar rotation for polar plot of color values
## used in seasonal concentration plots
prot=-180


## create 8-bit color table file for grass
cols=data.frame(id=0:104,col=apply(t(col2rgb(colR(105))),1,function(x) paste(x,collapse=":")))
write.table(cols,"data/out/grasscols.txt",col.names=F,row.names=F,quote=F)

## read in coast line
coast=readOGR("data/gshhs/","coast")


#######  Functions

## Long term summaries
seasconc <- function(x,return.Pc=T,return.thetat=T) {
  #################################################################################################
  ## Precipitation Concentration function
  ## This function calculates Precipitation Concentration based on Markham's (1970) technique as described in Schulze (1997)
  ## South Africa Atlas of Agrohydology and Climatology - R E Schulze, M Maharaj, S D Lynch, B J Howe, and B Melvile-Thomson
  ## Pages 37-38
  #################################################################################################
  ## x is a vector of precipitation quantities - the mean for each factor in "months" will be taken,
  ## so it does not matter if the data are daily or monthly, as long as the "months" factor correctly
  ## identifies them into 12 monthly bins, collapse indicates whether the data are already summarized as monthly means.
  #################################################################################################
  theta=seq(30,360,30)*(pi/180)                                       # set up angles for each month & convert to radians
  if(any(is.na(x))) { return(cbind(Pc=NA,thetat=NA))}
  if(return.Pc) {
    rt=sqrt(sum(x * cos(theta))^2 + sum(x * sin(theta))^2)    # the magnitude of the summation
    Pc=as.integer(round((rt/sum(x))*1000))}
  if(return.thetat){
    s1=sum(x*sin(theta),na.rm=T); s2=sum(x*cos(theta),na.rm=T)
    if(s1>=0 & s2>=0)  {thetat=abs((180/pi)*(atan(sum(x*sin(theta),na.rm=T)/sum(x*cos(theta),na.rm=T))))}
    if(s1>=0 & s2<=0)  {thetat=180-abs((180/pi)*(atan(sum(x*sin(theta),na.rm=T)/sum(x*cos(theta),na.rm=T))))}
    if(s1<=0 & s2<=0)  {thetat=180+abs((180/pi)*(atan(sum(x*sin(theta),na.rm=T)/sum(x*cos(theta),na.rm=T))))}
    if(s1<=0 & s2>=0)  {thetat=360-abs((180/pi)*(atan(sum(x*sin(theta),na.rm=T)/sum(x*cos(theta),na.rm=T))))}
    thetat=as.integer(round(thetat*10))
  }
  if(return.thetat&return.Pc) return(c(conc=Pc,theta=thetat))
  if(return.Pc)          return(Pc)
  if(return.thetat)  return(thetat)
}

### 
knitsDoc <- function(name) {
  library(knitr)
  knit(paste0(name, ".Rmd"), encoding = "utf-8")
  system(paste0("pandoc -o ", name, ".docx ", name, ".md --bibliography manuscript/biblio/MODCF.bib"))
}


#### Definition of regions to use for subsetting
regs=list(
  Cascades=extent(c(-122.8,-118,44.9,47)),
  Hawaii=extent(c(-156.5,-154,18.75,20.5)),
  Boliva=extent(c(-71,-63,-20,-15)),
  Venezuela=extent(c(-69,-59,0,7)),
  Venezuela2=extent(c(-72.5,-61,5,11)),
  Indonesia=extent(c(93.2,140,-14,12)),
  CFR=extent(c(17.5,29,-35,-29)),
  CFR2=extent(c(15,33,-35,-25)),
  Madagascar=extent(c(46,52,-17,-12))
)

fratio=function(ext) abs(ext@xmax-ext@xmin)/abs(ext@ymax-ext@ymin)

fratio(regs[["CFR2"]])
fratio(regs[["Indonesia"]])


## udpate paths because initGRASS gets it wrong for grass70...
ePATH <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(Sys.getenv("GISBASE"), "/bin:", 
                        Sys.getenv("GISBASE"), "/scripts", ifelse(nchar(ePATH) == 
                                                                    0, "", ":"), ePATH, sep = ""))
eLDPATH <- Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv(LD_LIBRARY_PATH = paste(Sys.getenv("GISBASE"), 
                                   "/lib:", ifelse(nchar(eLDPATH) == 0, "", ":"), 
                                   eLDPATH, sep = ""))
ePyPATH <- Sys.getenv("PYTHONPATH")
GrPyPATH <- paste(Sys.getenv("GISBASE"), "etc", "python", 
                  sep = "/")
Sys.setenv(PYTHONPATH = paste(GrPyPATH, ePyPATH, 
                              sep = ":"))

