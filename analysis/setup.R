
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