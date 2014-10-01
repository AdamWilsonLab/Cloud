# Load the libraries and set working directory
source("analysis/setup.R")

ncores=12
registerDoMC(ncores)

#sp=c("Lepidocolaptes","lacrymiger")
#fam="Furnariidae (Ovenbirds and Woodcreepers)"

## cock of the rock
sp=c("Rupicola","peruvianus")
fam="Cotingidae (Cotingas)"

sp2=paste0(sp,collapse="_")

### Load eBird data for 
## ebird data downloaded and unzipped from http://ebirddata.ornith.cornell.edu/downloads/erd/ebird_all_species/erd_western_hemisphere_data_grouped_by_year_v5.0.tar.gz
ebirddir="/mnt/data/jetzlab/Data/specdist/global/ebird/"

## select species for presences and 'absences'
taxon=read.csv(paste0(ebirddir,"/doc/taxonomy.csv"),na.strings="?")
taxon$TAXON_ORDER2=round(taxon$TAXON_ORDER)
## some duplicate species names marked with "/"
## Many SPECIES_NAME and GENUS_NAME are null but SCI_NAME is not
taxon$gensp=apply(do.call(rbind,lapply(strsplit(as.character(taxon$SCI_NAME),split="_|/"),function(x) x[1:2])),1,paste,collapse="_")
taxon=taxon%.%filter(FAMILY_NAME==fam)

## Select list (or single) taxon ids for use as 'presence'
sptaxon=taxon$TAXON_ORDER2[taxon$gensp%in%sp2]
## Select list of taxon ids for use as 'non-detection/absence'
nulltaxon=taxon$TAXON_ORDER2[taxon$TAXON_ORDER2!=sptaxon]

## create output folder
outputdir=paste0("output/sdm/",sp2,"/")
if(!file.exists(outputdir)) dir.create(outputdir,recursive=T)

## load region boundary
download.file(paste0("http://mol.cartodb.com/api/v2/sql?",
                  "q=SELECT%20ST_TRANSFORM(the_geom_webmercator,4326)%20as%20the_geom,%20seasonality%20FROM%20",
                  "get_tile('jetz','range','",
                  sp[1],"%20",sp[2],
                  "','jetz_maps')&format=shp&filename=expert"),
              destfile=paste0(outputdir,"expert.zip"))
unzip(paste0(outputdir,"expert.zip"),exdir=outputdir)


reg=readOGR(outputdir,"expert")
ereg=extent(reg)

## adjust bbox if desired
ereg@xmin=-81.4

## get species data
t1=system.time(spd<<-getebird(con=conn,sptaxon=sptaxon,nulltaxon=nulltaxon,region=reg))
coordinates(spd)=c("longitude","latitude")
projection(spd)="+proj=longlat +datum=WGS84 +ellps=WGS84"
spd@data[,c("lon","lat")]=coordinates(spd)  

## print a little summary
writeLines(paste(sp2," has ",nrow(spd),"rows"))

#########  Environmental Data
cf=stack(c("data/MCD09/MCD09_mean_01.tif","data/MCD09/MCD09_mean_07.tif","data/MCD09_deriv/intra.tif"))
names(cf)=c("CLDJAN","CLDJUL","CLDSEAS")
NAvalue(cf)=0
cf=crop(cf,ereg)
gain(cf)=.01

tvars=c(
  paste0("/mnt/data/jetzlab/Data/environ/global/worldclim/",
      c("prec_1","prec_7","bio_1","bio_12","bio_15","tmean_1","tmean_7","alt"),".bil"))
tenv=stack(crop(stack(tvars),ereg),cf)
## Covariate correlation
tcor=layerStats(tenv, "pearson",asSample=F, na.rm=T)
write.csv(tcor,file=paste0(outputdir,sp2,"_covariatecorrelation.csv"),row.names=F)


### Select environmental data
##  Create single stack
env=stack(tenv[[c("prec_1","prec_7","bio_15","bio_1","alt")]],cf[[c("CLDJAN","CLDJUL","CLDSEAS")]])
names(env)=c("PPTJAN","PPTJUL","PPTSEAS","MAT","ALT","CLDJAN","CLDJUL","CLDSEAS")

## rasterize species data to environmental grid
presences=rasterize(spd[spd$presence==1,],cf,fun="count",field="presence",background=0)
trials=rasterize(spd[spd$presence==0,],cf,fun="count",field="presence",background=0)
spr=stack(presences,trials)
names(spr)=c("presences","trials")

senv=scale(env)


data=cbind.data.frame(values(stack(senv,spr)),coordinates(senv),cell=cellFromXY(senv,xy=coordinates(senv)))
data=as.data.frame(data[!is.na(data[,"presences"]),])
data$trials=as.integer(data$trials)
data$presences=as.integer(data$presences)
data=na.omit(data)

## create 'fitting' dataset where there are observations
fdata=data[data$trials>0,]
fdata=data[data$trials>0&data$trials>=data$presences,]

## due to opportunistic observations, there are a few sites with more presences than trials, update those records so presences=trials
fdata$trials[fdata$presences>fdata$trials]=fdata$presences[fdata$presences>fdata$trials]

### Save model input data
save(env,senv,spr,data,fdata,file=paste0(outputdir,"modelinput.Rdata"))

#####################################################################################
#####################################################################################
load(paste0(outputdir,"modelinput.Rdata")))

nchains=3

mods=data.frame(
  model=c("m1","m2"),
  formula=c("~PPTJAN+PPTJUL+PPTSEAS+MAT",
            "~CLDJAN+CLDJUL+CLDSEAS+MAT"),
  name=c("Precipitation",
         "Cloud"))

## file to save output
fres=paste0(outputdir,"modeloutput.Rdata")

registerDoMC(ncores)

if(!file.exists(fres)){
res=foreach(m=1:nrow(mods)) %dopar% { 
  tres1=foreach(ch=1:nchains) %dopar% {
    tres2=hSDM.ZIB(
      suitability=as.character(mods$formula[m]),
      presences=fdata$presences,
      observability=~1,
      mugamma=0, Vgamma=1.0E6,
      gamma.start=0,
      trials=fdata$trials,
      data=fdata,
      burnin=10000, mcmc=10000, thin=10,
      beta.start=0,
      suitability.pred=data,
      mubeta=0, Vbeta=1.0E6,
      save.p=0,
      verbose=1,
      seed=round(runif(1,0,1000)))
}



