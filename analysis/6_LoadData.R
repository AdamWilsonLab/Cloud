## File to load data for use in SM.Rnw and 6_Figures


### Load data
print("Loading full raster layers")
cf_mean=raster("data/MCD09_deriv/meanannual.tif")
gain(cf_mean)=.01
cf_visseas=raster("data/MCD09_deriv/seas_visct.tif")
NAvalue(cf_visseas)=65534

inter=raster("data/MCD09_deriv/inter.tif")
gain(inter)=.01
intra=raster("data/MCD09_deriv/intra.tif")
gain(intra)=.01
NAvalue(inter)=0




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
