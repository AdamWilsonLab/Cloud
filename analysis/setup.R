
## libraries
library(rasterVis)
library(latticeExtra)
library(xtable)
library(texreg)
library(reshape2)
library(caTools)
library(rgeos)
library(raster)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(GGally)
library(scales)

library(magrittr)
library(tidyr)
library(knitr)
require(knitcitations)
## read in global coasts for nice plotting
library(maptools)
library(mapdata)

library(rgdal)
library(coda)
#library(spgrass6)
#install_github(repo="Rdatatable/data.table")
library(data.table)
library(AUC)
library(dismo)
library(redshift)
library(ncdf4)

library(dismo)
library(grid) # needed for arrow function


#library(devtools)

## install rasterAutocorr
library(devtools) 
#install_github("adammwilson/rasterAutocorr")
library(rasterAutocorr)

## update hSDM
#install_git("http://git.code.sf.net/p/hsdm/code", branch = "dev")
library(hSDM)

## register parallel backend
library(doMC)
ncores=12
registerDoMC(ncores)

rasterOptions(progress="text",maxmemory=1e6,tmpdir="data/tmp/",datatype="FLT4S")

## source all R files in R directory
sapply(list.files(pattern="[.]R$", path="R/", full.names=TRUE), source)

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

bgr=function(x,n=100,br=0,c1=c("blue","darkblue","grey"),c2=c("grey","darkred","red")){
  at=unique(c(seq(min(x,na.rm=T),max(x,na.rm=T),len=n)))
  #at=unique(c(seq(min(x,na.rm=T),br,len=floor(n/2)),seq(br,max(x,na.rm=T),len=ceiling(n/2))))
  bg=colorRampPalette(c1)
  gr=colorRampPalette(c2)
  return(list(at=at,col=c(bg(sum(at<br)),gr(sum(at>=br)))))
}

## set plotting parameters
my.theme = trellis.par.get()
my.theme$strip.background=list(col="transparent")
trellis.par.set(my.theme)


## Set polar rotation for polar plot of color values
## used in seasonal concentration plots
prot=-180


## create 8-bit color table file for grass
cols=data.frame(id=0:104,col=apply(t(col2rgb(colR(105))),1,function(x) paste(x,collapse=":")))
write.table(cols,"data/out/grasscols.txt",col.names=F,row.names=F,quote=F)

## read in coast line
coast=readOGR("data/gshhs/","coast")
fcoast=fortify(coast)

## coastline for ggplot
gcoast <- fortify(map_data('world'))#, plot = FALSE, fill = TRUE)
gcoast2 <- map_data('worldHires')#, plot = FALSE, fill = TRUE)


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
  Venezuela2=extent(c(-85,-55,0,12)),
  Venezuela3=extent(c(-85,-53,-7,11)),
  SouthAmerica=extent(c(-84,-33,-55,14)),
  Indonesia=extent(c(93.2,140,-14,12)),
  CFR=extent(c(17.5,29,-35,-29)),
  CFR2=extent(c(15,33,-35,-27.8)),
  Madagascar=extent(c(46,52,-17,-12)),
  Borneo=extent(c(107.5, 119.0918,-4.0835, 7.05)),
  EastAfrica=extent(c( 27,40.786056, -11.103410, 3.304579))

)

fratio=function(ext) abs(ext@xmax-ext@xmin)/abs(ext@ymax-ext@ymin)

abs(-180-180)/abs(84--60)
fratio(regs[["CFR2"]])
fratio(regs[["Venezuela2"]])


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

