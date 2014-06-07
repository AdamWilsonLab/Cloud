
## libraries
library(rasterVis)
library(latticeExtra)
library(xtable)
library(texreg)
library(reshape)
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
library(spgrass6)

## install rasterAutocorr
#library(devtools) 
#install_github("adammwilson/rasterAutocorr")
library(rasterAutocorr)

## register parallel backend
library(doMC)
registerDoMC(12)

rasterOptions(progress="text",maxmemory=1e9)


# temporary files will be written here:
datadir="/mnt/data2/projects/cloud/"

## color ramps
colR=colorRampPalette(c("#08306b","#0d57a1","#2878b8","#4997c9","#72b2d7","#a2cbe2","#c7dcef","#deebf7","#f7fbff"))
bgr=colorRampPalette(c("#0000ff","#00ff00","#ff0000"))
bgyrp=colorRampPalette(c("blue","darkgreen","goldenrod","red","purple"))


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
  CFR=extent(c(17.75,22.5,-34.8,-32.6)),
  Madagascar=extent(c(46,52,-17,-12))
)




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

