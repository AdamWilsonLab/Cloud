### This script performs additional pre-processing of datasets used in this project

source("analysis/setup.R")




## Global coastline
land=readShapePoly("/mnt/data/jetzlab/Data/environ/global/gshhg/GSHHS_shp/c/GSHHS_c_L1.shp",force_ring=TRUE)
projection(land)="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
CP <- as(extent(-180, 180, -60, 84), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(land))
coast=as(land[land$area>50,],"SpatialLines")
land <- gIntersection(land, CP, byid=F)
coast <- gIntersection(coast, CP, byid=F)

hland=readShapePoly("/mnt/data/jetzlab/Data/environ/global/gshhg/GSHHS_shp/i/GSHHS_i_L1.shp",force_ring=TRUE)
projection(hland)="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#hland=hland[hland$area>50,]
hcoast=as(hland,"SpatialLinesDataFrame")
hcoast <- gIntersection(hcoast, CP, byid=F)
hcoast<-as(hcoast,"SpatialLinesDataFrame")
hcoast@data=data.frame(ID=1)
writeOGR(hcoast,"data/gshhs/","coast",driver="ESRI Shapefile",overwrite=T)

