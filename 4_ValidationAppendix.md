Appendix A: Methods & Validation
=======================

This appendix includes details of the validation procedure.



```r
#,setup,echo=F,results='hide',message=F}
## some setup

uploadimages=F  # upload all images to imgur.com for easy viewing on github
opts_knit$set(progress = TRUE, verbose = TRUE,root.dir="../../",cache=!uploadimages)
#if(uploadimages) opts_knit$set(upload.fun =imgur_upload)
#opts_chunk$set(fig.width=12, fig.height=8, cache=!uploadimages)
getwd()
```

```
## [1] "/media/data/Cloud/manuscript/appendices"
```

```r
## libraries
library(rasterVis)
```

```
## Loading required package: raster
## Loading required package: sp
## Loading required package: lattice
## Loading required package: latticeExtra
## Loading required package: RColorBrewer
## Loading required package: hexbin
## Loading required package: grid
```

```r
library(latticeExtra)
library(xtable)
library(texreg)
```

```
## Version:  1.31
## Date:     2014-02-08
## Author:   Philip Leifeld (University of Konstanz)
## 
## Attaching package: 'texreg'
## 
## The following object is masked from 'package:raster':
## 
##     extract
```

```r
library(reshape)
```

```
## Loading required package: plyr
## 
## Attaching package: 'plyr'
## 
## The following object is masked from 'package:raster':
## 
##     count
## 
## 
## Attaching package: 'reshape'
## 
## The following objects are masked from 'package:plyr':
## 
##     rename, round_any
## 
## The following object is masked from 'package:raster':
## 
##     expand
```

```r
library(caTools)
library(rgeos)
```

```
## rgeos version: 0.2-19, (SVN revision 394)
##  GEOS runtime version: 3.3.8-CAPI-1.7.8 
##  Polygon checking: TRUE
```

```r
library(raster)
library(plyr)
library(knitr)
require(knitcitations)
```

```
## Loading required package: knitcitations
## Loading required package: bibtex
## 
## Attaching package: 'knitcitations'
## 
## The following object is masked from 'package:utils':
## 
##     cite
```

```r
## read in global coasts for nice plotting
library(maptools)
```

```
## Loading required package: foreign
## Checking rgeos availability: TRUE
## 
## Attaching package: 'maptools'
## 
## The following object is masked from 'package:xtable':
## 
##     label
```

```r
library(rgdal)
```

```
## rgdal: version: 0.8-10, (SVN revision 478)
## Geospatial Data Abstraction Library extensions to R successfully loaded
## Loaded GDAL runtime: GDAL 1.10.0, released 2013/04/24
## Path to GDAL shared files: /usr/share/gdal/1.10
## Loaded PROJ.4 runtime: Rel. 4.8.0, 6 March 2012, [PJ_VERSION: 480]
## Path to PROJ.4 shared files: (autodetected)
```


```
## Warning: cannot open file '../../data/validation/cldm.csv': No such file
## or directory
```

```
## Error: cannot open the connection
```

```
## Error: Cannot open file
```

```
## Error: object 'cldm' not found
```

```
## Error: object 'cldm' not found
```

```
## Error: object 'cldm' not found
```

```
## Error: object 'cldm' not found
```

```
## Error: object 'cldm' not found
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function '%in%': Error: object 'cldm' not found
```

```
## Error: object 'cldm' not found
```

```
## Error: object 'cldm' not found
```

```
## Error: object 'cldml' not found
```


```
## Error: error in evaluating the argument 'x' in selecting a method for function 'unique': Error: object 'cldm' not found
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'unique': Error: object 'cldm' not found
```

















