
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

## register parallel backend
library(doMC)
registerDoMC(12)

## Working Directory
setwd("/media/data/Cloud")


## Load coast object for easy plotting
land=readShapePoly("/mnt/data/jetzlab/Data/environ/global/gshhg/GSHHS_shp/c/GSHHS_c_L1.shp",force_ring=TRUE)
projection(land)="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
CP <- as(extent(-180, 180, -60, 84), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(land))
coast=as(land[land$area>50,],"SpatialLines")
land <- gIntersection(land, CP, byid=F)
coast <- gIntersection(coast, CP, byid=F)


#######  Functions

## Long term summaries
seasconc <- function(x,return.Pc=T,return.thetat=F) {
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
  if(sum(is.na(x))==12) { return(cbind(Pc=NA,thetat=NA)) ; stop}
  if(return.Pc) {
    rt=sqrt(sum(x * cos(theta))^2 + sum(x * sin(theta))^2)    # the magnitude of the summation
    Pc=as.integer(round((rt/sum(x))*100))}
  if(return.thetat){
    s1=sum(x*sin(theta),na.rm=T); s2=sum(x*cos(theta),na.rm=T)
    if(s1>=0 & s2>=0)  {thetat=abs((180/pi)*(atan(sum(x*sin(theta),na.rm=T)/sum(x*cos(theta),na.rm=T))))}
    if(s1>0 & s2<0)  {thetat=180-abs((180/pi)*(atan(sum(x*sin(theta),na.rm=T)/sum(x*cos(theta),na.rm=T))))}
    if(s1<0 & s2<0)  {thetat=180+abs((180/pi)*(atan(sum(x*sin(theta),na.rm=T)/sum(x*cos(theta),na.rm=T))))}
    if(s1<0 & s2>0)  {thetat=360-abs((180/pi)*(atan(sum(x*sin(theta),na.rm=T)/sum(x*cos(theta),na.rm=T))))}
    thetat=as.integer(round(thetat))
  }
  if(return.thetat&return.Pc) return(c(conc=Pc,theta=thetat))
  if(return.Pc)          return(Pc)
  if(return.thetat)  return(thetat)
}

