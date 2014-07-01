## File to load data for use in SM.Rnw and 6_Figures


### Load data
print("Loading full raster layers")
cf_mean=readAll(raster("data/MCD09_deriv/MCD09_meanannual.tif"))
cf_visseas=readAll(raster("data/MCD09_deriv/seas_visct.tif"))
inter=readAll(raster("data/MCD09_deriv/inter.tif"))
intra=readAll(raster("data/MCD09_deriv/intra.tif"))
NAvalue(inter)=0

### Load validation data
cldm=read.csv("data/validation/cldm.csv",row.names=NULL)
st=readOGR("data/validation","stations")
st$era=factor(st$era,levels=c("Pre-MODIS","Full"),labels=c("Pre-MODIS","Full"),ordered=T)

## month factors
cldm$month2=factor(cldm$month,labels=month.name,ordered=T)

### Drop valitation station-months with fewer than 20 years of data for full record or less than 10 years for MODIS-era record
cldm$cld_all[cldm$cldn_all<20]=NA
cldm$cldsd_all[cldm$cldn_all<20]=NA

cldm$cld[cldm$cldn<10]=NA
cldm$cldsd[cldm$cldn<10]=NA

cldm$seas=factor(cldm$seas,levels=c("DJF","MAM","JJA","SON"),ordered=T)

## compute seasonal means
cldml=melt(cldm,id.vars=c("StaID","lat","lon","lulcc","biome","seas"),measure.vars=c("cld_all","cldsd_all","cldn_all","cld","cldsd","cldn","MCD09_meanb16","MCD09_meanb5","MCD09_sdb16"))
clds=dcast(cldml,StaID+lat+lon+lulcc+biome+seas~variable,value.var="value",fun=mean,na.rm=T)

# get residuals of simple linear model
cldm$resid=NA
cldm$resid[!is.na(cldm$cld_all)&!is.na(cldm$MCD09_meanb16)]=residuals(lm(MCD09_meanb16~cld_all,data=cldm[!is.na(cldm$cld_all)&!is.na(cldm$MCD09_meanb16),]))

# get residuals of simple linear model
cldm$difm_all=cldm$MCD09_meanb16-cldm$cld_all
cldm$difm=cldm$MCD09_meanb16-cldm$cld



########################################
### Some stats 
## number of stations retained
nstation_all=length(unique(cldm$StaID[!is.na(cldm$cld_all)]))
nstation_mod=length(unique(cldm$StaID[!is.na(cldm$cld)]))

# approximate size of M*D09GA archive - get total size for one day from the USGS website
size=as.numeric(
  sub("M","",
      grep("[0-9]*M$",
           scan("http://e4ftl01.cr.usgs.gov/MOLT/MOD09GA.005/2000.04.30/",what="char"),
           value=T)))
## extract all filesizes in MB (all the HDFs) and sum them and covert to TB for the length of the full record
if(!file.exists("data/out/SizeOfMOD09GAArchive.Rdata")){
  tsize=sum(size)/1024/1024*as.numeric(as.Date("2014-03-31")-as.Date("2000-02-24"))*2
  save(tsize,file="data/out/SizeOfMOD09GAArchive.Rdata")
}
load("data/out/SizeOfMOD09GAArchive.Rdata")

## get HPD intervals
#hpd_resid=HPDinterval(as.mcmc(cldm$resid))#,c(0.025,0.25,0.5,0.75,0.975),na.rm=T)
#hpd_cld=HPDinterval(as.mcmc(cldm$difm))#,c(0.025,0.25,0.5,0.75,0.975),na.rm=T)
#hpd_cld_all=HPDinterval(as.mcmc(cldm$difm_all))#,c(0.025,0.25,0.5,0.75,0.975),na.rm=T)
