# Load the libraries and set working directory
source("analysis/setup.R")

ncores=12
registerDoMC(ncores)

sp=c("Lepidocolaptes","lacrymiger")
fam="Furnariidae (Ovenbirds and Woodcreepers)"

## cock of the rock
#sp=c("Rupicola","peruvianus")
#fam="Cotingidae (Cotingas)"

sp2=paste0(sp,collapse="_")

## create output folder
outputdir=paste0("output/sdm/",sp2,"/")
if(!file.exists(outputdir)) dir.create(outputdir,recursive=T)


## file to hold model input data
fmodelinput=paste0(outputdir,"/",sp2,"_modelinput.nc")

##########################################################

### Load eBird data for 
## ebird data downloaded and unzipped from http://ebirddata.ornith.cornell.edu/downloads/erd/ebird_all_species/erd_western_hemisphere_data_grouped_by_year_v5.0.tar.gz
ebirddir="/mnt/data/jetzlab/Data/specdist/global/ebird/"

## select species for presences and 'absences'
taxon=read.csv(paste0(ebirddir,"/doc/taxonomy.csv"),na.strings="?",stringsAsFactors=F)
taxon$TAXON_ORDER2=floor(taxon$TAXON_ORDER)
## some duplicate species names marked with "/"
## Many SPECIES_NAME and GENUS_NAME are null but SCI_NAME is not, merge them (dropping subspecies,etc)
## take first word of SCI_NAME (before "_") and call that genus
taxon$GENUS_NAME[is.na(taxon$GENUS_NAME)]=do.call(rbind,lapply(
  strsplit(as.character(taxon$SCI_NAME[is.na(taxon$GENUS_NAME)]),split="_|/"),function(x) x[1]))
## take second word of SCI_NAME (before "_") and call that species
taxon$SPECIES_NAME[is.na(taxon$SPECIES_NAME)]=do.call(rbind,lapply(
  strsplit(as.character(taxon$SCI_NAME[is.na(taxon$SPECIES_NAME)]),split="_|/"),function(x) x[2]))
## put them together
taxon$gensp=paste(taxon$GENUS_NAME,taxon$SPECIES_NAME,sep="_")

# Filter by family?
#taxon=taxon%.%filter(FAMILY_NAME==fam)

## Select list (or single) taxon ids for use as 'presence'
taxon[taxon$gensp%in%sp2,]
sptaxon=taxon$TAXON_ORDER2[taxon$gensp%in%sp2]

## Select list of taxon ids for use as 'non-detection/absence' - if desired
nulltaxon=taxon$TAXON_ORDER2[taxon$TAXON_ORDER2!=sptaxon]


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
t1=system.time(spd_all<<-getebird(con=conn,sptaxon=sptaxon,nulltaxon=NULL,region=reg))
writeLines(paste("eBird extraction for",sp2," took ",round(t1[[3]]/60,2),"seconds and returned",nrow(spd_all),"records"))

## trim observations by observer effort, distance travelled, etc.
cdur=4*60
cdis=5
care=500

ggplot(spd_all,aes(y=duration_minutes,x=effort_distance_km,colour=as.factor(presence),order=as.factor(presence)))+
  geom_point()+scale_x_log10()+
  geom_vline(xintercept=cdis)+geom_hline(yintercept=cdur)

#spd_woodcreeper=spd
spd=filter(spd_all,duration_minutes<=cdur&(effort_distance_km<=cdis|effort_area_ha<=care))

table(presence=spd_all$presence>0,distance=cut(spd_all$effort_distance_km,c(0,1,5,10,50,100,Inf)))
table(presence=spd_all$presence>0,time=cut(spd_all$duration_minutes,c(0,30,60,120,240,480,1440)))
table(durationNA=is.na(spd_all$duration_minutes),distanceNA=is.na(spd_all$effort_distance_km),areaNA=is.na(spd_all$effort_area_ha))

table(presence=spd_all$presence,filter=spd_all$duration_minutes>cdur&(spd_all$effort_distance_km>cdis|spd_all$effort_area_ha>care))


coordinates(spd)=c("longitude","latitude")
projection(spd)="+proj=longlat +datum=WGS84 +ellps=WGS84"
spd@data[,c("lon","lat")]=coordinates(spd)  


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

## define metadata that will be embeded in output object
meta=list(institution="Map of Life, Yale University, New Haven, CT",
          source="Species Distributions",
          comment="Adam M. Wilson (adam.wilson@yale.edu)",
          species=sp2)

hSDM.ncWriteInput(env,points=spd,ncfile=fmodelinput,overwrite=T,meta=meta)

##############################################################
## import data
data=hSDM.ncReadInput(fmodelinput)
## select data for fitting
## due to opportunistic observations, there are a few sites with more presences than trials,
data$fit=ifelse(data$trials>0&data$trials>=data$presences,1,0)
## create 'fitting' dataset where there are observations
fdata=data[data$fit>0,]


nchains=3

mods=data.frame(
  model=c("m1","m2"),
  formula=c("~PPTJAN+PPTJUL+PPTSEAS+MAT",
            "~CLDJAN+CLDJUL+CLDSEAS+MAT"),
  name=c("Precipitation",
         "Cloud"))

registerDoMC(ncores)

burnin=10000
mcmc=10000
thin=10

foreach(m=1:nrow(mods)) %dopar% { 
  tm=P.hSDM.ZIB(
    suitability=as.character(mods$formula[m]),
    presences=fdata$presences,
    observability=~1,
    mugamma=0, Vgamma=1.0E6,
    gamma.start=0,
    trials=fdata$trials,
    data=fdata,
    burnin=burnin, mcmc=mcmc, thin=thin,
    beta.start=0,
    suitability.pred=data,
    mubeta=0, Vbeta=1.0E6,
    save.p=0,
    verbose=1,
    seed=round(runif(1,0,1e6)))
  
  meta=list(modelname=as.character(mods$name[m]),species=sp2)
  outfile=paste0(outputdir,"/",sp2,"_modeloutput_",mods$name[m],".nc")
  source("/media/data/repos/hsdm-code/R/hSDM.ncWriteOutput.R")  
    hSDM.ncWriteOutput(results=tm,file=outfile,meta=meta,verbose=T,autocor=T)    
}
 
