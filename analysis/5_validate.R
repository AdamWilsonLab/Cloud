### Script to download and process the NDP-026D station cloud dataset
### to validate MODIS cloud frequencies

source("analysis/setup.R")

## Data available here http://cdiac.ornl.gov/epubs/ndp/ndp026d/ndp026d.html

download=F  #download data?
## Get station locations
if(download)   system("wget -N -nd http://cdiac.ornl.gov/ftp/ndp026d/cat01/01_STID -P data/validation/NDP026D/")

st=read.table("data/validation/NDP026D/01_STID",skip=1)
colnames(st)=c("StaID","LAT","LON","ELEV","ny1","fy1","ly1","ny7","fy7","ly7","SDC","b5c")
st$lat=st$LAT/100
st$lon=st$LON/100
st$lon[st$lon>180]=st$lon[st$lon>180]-360
st=st[,c("StaID","ELEV","lat","lon")]
colnames(st)=c("id","elev","lat","lon")
write.csv(st,"stations.csv",row.names=F)
coordinates(st)=c("lon","lat")
projection(st)="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
st@data[,c("lon","lat")]=coordinates(st)

## download data
if(download){
    system("wget -N -nd ftp://cdiac.ornl.gov/pub/ndp026d/cat67_78/* -A '.tc.Z' -P data/validation/NDP026D/")
    system("gunzip data/*.Z")
}

## define FWF widths
f162=c(5,5,4,7,7,7,4) #format 162
c162=c("StaID","YR","Nobs","Amt","Fq","AWP","NC")

## use monthly timeseries
cld=do.call(rbind.data.frame,mclapply(sprintf("%02d",1:12),function(m) {
  d=read.fwf(list.files("data/validation/NDP026D/",pattern=paste("MNYDC.",m,".tc$",sep=""),full=T),skip=1,widths=f162)
  colnames(d)=c162
  d$month=as.numeric(m)
  print(m)
  return(d)}
))

## add lat/lon
cld[,c("lat","lon")]=coordinates(st)[match(cld$StaID,st$id),]

## drop missing values
cld=cld[,!grepl("Fq|AWP|NC",colnames(cld))]
cld$Amt[cld$Amt<0]=NA
cld$Amt=cld$Amt/100
cld=cld[!is.na(cld$Amt),]

## table of stations with > 20 observations per month
#dcast(cld,StaID~YR,value.var="Nobs",fun=sum)
mtab=cld %.% group_by(StaID,month) %.% summarise(count=sum(Nobs>20,na.rm=T)) 
mtab2=mtab %.% group_by(StaID) %.% summarise(count=any(count>10))
stem(mtab$count)

## calculate means and sds for full record (1970-2009)
stem(cld$Nobs)
Nobsthresh=20 #minimum number of observations to include 

cldm=do.call(rbind.data.frame,by(cld,list(month=as.factor(cld$month),StaID=as.factor(cld$StaID)),function(x){
  data.frame(
    month=x$month[1],
    StaID=x$StaID[1],
    cld_all=mean(x$Amt[x$Nobs>=Nobsthresh],na.rm=T),  # full record
    cldsd_all=sd(x$Amt[x$Nobs>=Nobsthresh],na.rm=T),
    cldn_all=length(x$Amt[x$Nobs>=Nobsthresh]),
    cld=mean(x$Amt[x$YR>=2000&x$Nobs>=Nobsthresh],na.rm=T), #only MODIS epoch
    cldsd=sd(x$Amt[x$YR>=2000&x$Nobs>=Nobsthresh],na.rm=T),
    cldn=length(x$Amt[x$YR>=2000&x$Nobs>=Nobsthresh]))}))

cldm[,c("lat","lon")]=coordinates(st)[match(cldm$StaID,st$id),c("lat","lon")]



## add the EarthEnvCloud data to cld
mod09_mean=stack((list.files("data/MCD09/",pattern="MCD09_mean_[0-9]*[.]tif",full=T)))
NAvalue(mod09_mean)=65535
gain(mod09_mean)=.01
names(mod09_mean)=month.name

mod09_sd=stack(list.files("data/MCD09/",pattern="MCD09_sd_[0-9]*[.]tif",full=T))
NAvalue(mod09_sd)=65535
gain(mod09_sd)=.01
names(mod09_sd)=month.name


## overlay the data with 32km diameter (16km radius) buffer
## buffer size from Dybbroe, et al. (2005) doi:10.1175/JAM-2189.1.
buf=16000
bins=cut(st$lat,quantile(st$lat,seq(0,1,len=20)))
## should cut by quantile instead to produce more even jobs...

rerun=F
if(rerun&file.exists("valid.csv")) file.remove(list.files("data/validation/",pattern="valid_tiled",full=T))

foreach(lb=levels(bins)) %dopar% {
  l=which(bins==lb)
  ## mean
  td1=raster::extract(mod09_mean,st[l,],buffer=buf,fun=mean,na.rm=T,df=T)
  td1$id=st$id[l]
  td1$type="MCD09_mean"
  ## std
  td2=raster::extract(mod09_sd,st[l,],buffer=buf,fun=mean,na.rm=T,df=T)
  td2$id=st$id[l]
  td2$type="MCD09_sd"
  print(lb)
  write.csv(rbind(td1,td2),paste0("data/validation/valid_tiled_",gsub("\\(|\\]|,","",lb),".csv"),col.names=T,quote=F,row.names=F)
}


## read it back in
mod09st=do.call(rbind.data.frame,lapply(list.files("data/validation/",pattern="valid_tiled",full=T),function(f)
  read.csv(f)[,-c(1)]))
mod09stl=melt(mod09st,id.vars=c("id","type"))
colnames(mod09stl)[grep("variable",colnames(mod09stl))]="month"
mod09stl$value[mod09stl$value<0]=NA
mod09stl=dcast(mod09stl,id+month~type,value="value")

## add it to cld
cldm$monthname=month.name[cldm$month]
cldm$MCD09_mean=mod09stl$MCD09_mean[match(paste(cldm$StaID,cldm$monthname),paste(mod09stl$id,mod09stl$month))]
cldm$MCD09_sd=mod09stl$MCD09_sd[match(paste(cldm$StaID,cldm$monthname),paste(mod09stl$id,mod09stl$month))]


## LULC
#system(paste("gdalwarp -r near -co \"COMPRESS=LZW\" -tr ",paste(res(mod09),collapse=" ",sep=""),
#             "-tap -multi -t_srs \"",   projection(mod09),"\" /mnt/data/jetzlab/Data/environ/global/landcover/MODIS/MCD12Q1_IGBP_2005_v51.tif ../modis/mod12/MCD12Q1_IGBP_2005_v51.tif"))
lulc=raster("/mnt/data/personal/adamw/projects/interp/data/modis/mod12/MCD12Q1_IGBP_2005_v51.tif")
require(plotKML); data(worldgrids_pal)  #load IGBP palette
IGBP=data.frame(ID=0:16,col=worldgrids_pal$IGBP[-c(18,19)],stringsAsFactors=F)
IGBP$class=rownames(IGBP);rownames(IGBP)=1:nrow(IGBP)
levels(lulc)=list(IGBP)
## function to get modal lulc value
Mode <- function(x) {
      ux <- na.omit(unique(x))
        ux[which.max(tabulate(match(x, ux)))]
      }

lulcst=raster::extract(lulc,st,fun=Mode,buffer=buf,df=T)
colnames(lulcst)=c("id","lulc")
lulcst$StaID=st$id
## add it to cld
cldm$lulc=lulcst$lulc[match(cldm$StaID,lulcst$StaID)]
cldm$lulcc=IGBP$class[match(cldm$lulc,IGBP$ID)]


### Add biome data
biome=readOGR("data/src/teow/","biomes")
projection(biome)=projection(st)
#st$biome=over(st,biome,returnList=F)$BiomeID
dists=apply(gDistance(st,biome,byid=T),2,which.min)
st$biomec=biome$code[dists]
st$realm=biome$realm[dists]
st$biome=biome$biome[dists]

cldm$biomec=st$biomec[match(cldm$StaID,st$id)]
cldm$realm=st$relam[match(cldm$StaID,st$id)]
cldm$biome=st$biome[match(cldm$StaID,st$id)]

cldm$seas=ifelse(cldm$month%in%c(12,1,2),"DJF",
                 ifelse(cldm$month%in%3:4,"MAM",
                        ifelse(cldm$month%in%5:8,"JJA",
                               ifelse(cldm$month%in%9:11,"SON",NA))))

## add indicator for MODIS era to st
mtab=do.call(rbind.data.frame,by(cldm,cldm$StaID,function(x) c(All=any(!is.na(x$cld_all)),MODIS=any(!is.na(x$cld)))))
mtab$era=ifelse(mtab[,1]&mtab[,2],"Full",
                ifelse(mtab[,1]&!mtab[,2],"Pre-MODIS",
                       ifelse(!mtab[,1]&mtab[,2],"MODIS","None")))
mtab$id=rownames(mtab$id)
st$era=as.factor(mtab$era[match(as.character(st$id),rownames(mtab))])


## write out the tables
write.csv(cld,file="data/validation/cld.csv",row.names=F)
write.csv(cldm,file="data/validation/cldm.csv",row.names=F)
writeOGR(st,dsn="data/validation/",layer="stations",driver="ESRI Shapefile",overwrite_layer=T)


#########################################################################
cldm=read.csv("data/validation/cldm.csv")
st=readOGR("data/validation","stations")
st$era=factor(st$era,levels=c("Pre-MODIS","Full"),labels=c("Pre-MODIS","Full"),ordered=T)

## set factor ordering
cldm$monthname=factor(cldm$monthname,ordered=T,levels=month.name)
cldm$seas=factor(cldm$seas,ordered=T,levels=c("DJF","MAM","JJA","SON"))

## add elevation to cldm
cldm$elev=st$elev[match(cldm$StaID,st$id)]

## subset cldm to include only stations with full 10 years of data for 2000-2009 and >= 20 years for 1970-2009
cldm[cldm$cldn<10,c("cld","cldsd")]=NA
cldm[cldm$cldn_all<20,c("cld_all","cldsd_all")]=NA

## compute seasonal means
cldml=melt(cldm,id.vars=c("StaID","lat","lon","lulcc","biome","seas"),measure.vars=c("cld_all","cldsd_all","cldn_all","cld","cldsd","cldn","MCD09_mean","MCD09_sd"))
clds=dcast(cldml,StaID+lat+lon+lulcc+biome+seas~variable,value.var="value",fun=mean,na.rm=T)

# get residuals of simple linear model
cldm$resid=NA
cldm$resid[!is.na(cldm$cld_all)&!is.na(cldm$MCD09_mean)]=residuals(lm(MCD09_mean~cld_all,data=cldm[!is.na(cldm$cld_all)&!is.na(cldm$MCD09_mean),]))

# get residuals of simple linear model
cldm$difm_all=cldm$MCD09_mean-cldm$cld_all
cldm$difm=cldm$MCD09_mean-cldm$cld


source("R/lm_summary.R")
t1=cldm %.% group_by(monthname) %.% summarise(
  Mean=mean(MCD09_mean,na.rm=T),
  n=lm_summary(cld,MCD09_mean,"n"),
  R2=lm_summary(cld,MCD09_mean,"rs"),
  RMSE=lm_summary(cld,MCD09_mean,"rmse"))
t2=cldm %.% group_by(seas) %.% summarise(
  Mean=mean(MCD09_mean,na.rm=T),
  n=lm_summary(cld,MCD09_mean,"n"),
  R2=lm_summary(cld,MCD09_mean,"rs"),
  RMSE=lm_summary(cld,MCD09_mean,"rmse"))
t3=data.frame(
  "Annual",
  Mean=mean(cldm$MCD09_mean,na.rm=T),
  n=lm_summary(cldm$cld,cldm$MCD09_mean,"n"),
  R2=lm_summary(cldm$cld,cldm$MCD09_mean,"rs"),
  RMSE=lm_summary(cldm$cld,cldm$MCD09_mean,"rmse"))
colnames(t1)[1]="Month/Season"
colnames(t2)[1]="Month/Season"
colnames(t3)[1]="Month/Season"

## combine to a single table
vtable=rbind.data.frame(t3,t2,t1)

print(xtable(vtable,digits=2,#caption="Summary of validation data by month and season"),
  "html",format.args=list(big.mark = ",", decimal.mark = "."),include.rownames=F,file="manuscript/validtable.html"))


cldm %>% group_by(StaID) %>% summarise(
  MCD09=mean(MCD09_mean,na.rm=T),
  cld=mean(cld,na.rm=T),
  elev=mean(elev,na.rm=T)) %>%
  summarize(
    cld_R2=lm_summary(cld,elev,"rs"),
    MCD09_R2=lm_summary(MCD09,elev,"rs"),
    cld_rho=cor.test(cld,elev,method="spearman",alternative="two.sided",continuity=T)$estimate,
    cld_n=nrow(na.omit(cbind(cld,elev))),
    MCD09_rho=cor.test(MCD09,elev,method="spearman",alternative="two.sided",continuity=T)$estimate,
    MCD09_n=nrow(na.omit(cbind(MCD09,elev))))


########################
### Temporal stability
### Compare two time periods
lm_all1=lm(cld_all~MCD09_mean,data=cldm)
lm_mod=lm(cld~MCD09_mean,data=cldm)

mods=list("1970-2009"=lm_all1,"2000-2009"=lm_mod)

screenreg(mods,digits=2,single.row=T,custom.model.names=names(mods),custom.coef.names = c("Intercept", "MODCF"))

texreg(mods,custom.model.names = names(mods),custom.coef.names = c("Intercept", "MODCF"),
       single.row = T,caption="Comparison of validation models for full station record (1970-2009) and MODIS era (2000-2009). Stations were included if they had at least 20 years of data for full record or 10 years for MODIS-era record")#,caption.placement='top')




### Spatial plot
ggplot(cldm,aes(x=lon,y=lat,colour=difm_all,order=desc(abs(difm_all))))+
  scale_colour_gradient2(low="blue",mid="grey",high="red")+
  facet_wrap(~monthname)+
  geom_point()+
  coord_equal()

ggplot(cldm,aes(x=difm,y=lat,col=difm))+
  scale_colour_gradient2(low="blue",mid="grey",high="red",na.value="transparent",name=expression(paste(Delta,"CLD")))+
  facet_grid(month~.)+
  ylab("")+
  xlab(expression(paste(Delta,"CLD")))+
  geom_point(size=.75)+
  geom_vline(xintercept=0)+
  stat_summary(geom = "smooth")


###############
## Calculate seasonal differences
clds$difs=clds$MCD09_mean-clds$cld_all
clds$resid=NA
clds$resid[!is.na(clds$cld_all)&!is.na(clds$MCD09_mean)]=
  residuals(lm(MCD09_mean~cld_all,data=clds[!is.na(clds$cld_all)&!is.na(clds$MCD09_mean),]))


## Add panel ID
clds$panelid1=factor(clds$seas, labels=letters[1:4])
clds$panelid2=factor(clds$seas, labels=letters[5:8])


gp=gpar(fontsize=14, col="black",fontface="bold")

pmap=ggplot(clds,aes(x=lon,y=lat,colour=resid,order=desc(abs(resid))))+
  scale_colour_gradient2(low="blue",mid="grey",high="red",na.value="transparent",name="CF Residuals",guide=F)+
  facet_grid(seas~.)+
  geom_point(size=.75)+
  geom_path(data=gcoast,mapping=aes(x=long,y=lat,group=group,order=order),col="black",size=.2)+
  ylab("Latitude")+xlab("Longitude")+
  geom_text(aes(x=-180,y=80,label=panelid1),col="black",stat="unique",gp=gp)+
  coord_equal()+
  theme(strip.text.x = element_blank(),strip.text.y = element_blank(),strip.background = element_blank())

pzonal=ggplot(clds,aes(x=resid,y=lat,col=resid))+
  scale_colour_gradient2(low="blue",mid="grey",high="red",na.value="transparent",name="CF Residuals")+
  facet_grid(seas~.)+
  ylab("")+
  xlab("CF Residuals")+
  geom_point(size=.75)+
  geom_vline(xintercept=0)+
  geom_text(aes(x=-58,y=72,label=panelid2) ,col="black" ,stat="unique",gp=gp)+
  stat_summary(geom = "smooth")

png(width=2400,height=1600,res=300,pointsize=16,
    type="cairo-png",file="manuscript/figures/validationMap.png")
print(pmap, vp = viewport(width = .7, height = 1, x=.35,y=.5))
print(pzonal, vp = viewport(width = .5, height = 1, x = .75, y = 0.5))
dev.off()


#################
#####  LULC plots
ggplot(clds,aes(x=resid,group=seas,col=seas))+
  scale_colour_discrete(name="Season")+
  facet_wrap(~lulcc)+
  ylab("")+
  xlab("Cloud Frequency Residuals")+
  geom_density(size=.75)+
  geom_vline(xintercept=0)

plulc=
  ggplot(clds,aes(y=resid,x=seas,col=seas,fill=seas))+
  scale_colour_discrete(guide=F)+
  scale_fill_discrete(guide=F)+  
  facet_wrap(~lulcc)+
  xlab("Season")+
  ylab("Cloud Frequency Residuals")+
  geom_jitter(size=1,col=grey(.6))+
  geom_violin(size=.75,scale="count",alpha=.5)+
  geom_hline(yintercept=0,linetype="dashed",col=grey(.2))


png(width=3000,height=1600,res=300,pointsize=16,
    type="cairo-png",file="manuscript/figures/validationLULC.png")
plulc
dev.off()


ggplot(clds[clds$lulcc=="Snow and ice ",],aes(x=cld_all,y=MCD09_mean))+
  facet_grid(lulcc~seas)+
  geom_point(size=.75)+
  geom_abline()+
  stat_smooth(method="lm",col="red")+
  coord_equal()+
  theme(strip.text.y = element_text(angle=0))


## Table of LULC residuals

cldm %>% group_by(lulcc) %>% summarise(
  resid=mean(resid,na.rm=T),
  residsd=sd(resid,na.rm=T))

## Write out seasonal validation data
head(clds)
write.csv(select(clds,-c(panelid1,panelid2)),"output/SeasonalValidation.csv",row.names=F)


##############
### Extract regional transects
cids=c(10866,10980,11130,11135,11138,11146,11210,11212,13014)

cids=c(10858,10865,10866,10980,11130,11135,11138,11204,11207,16040) #11138,11146,11210,11212,13014)
ts=cldm[cldm$StaID%in%cids,]
ts$trans=match(ts$StaID,cids)
ts$cld_allpsd=ts$cld_all+ts$cldsd_all
ts$cld_allmsd=ts$cld_all-ts$cldsd_all
ts$mod09msd=ts$MCD09_mean-ts$MCD09_sd
ts$mod09psd=ts$MCD09_mean+ts$MCD09_sd
ts$mod09msd=ts$MCD09_mean-ts$MCD09_sd
ts$mod09psd=ts$MCD09_mean+ts$MCD09_sd

tsl=melt(ts,id.vars=c("month","trans"),measure.vars=c("cld_all","cld_allmsd","cld_allpsd","MCD09_mean","MCD09_sd"))

#as("SpatialLinesDataFrame")
xyplot(value~trans|month,groups=variable,data=tsl,type="l",auto.key=T,ylim=c(0,100))

xyplot(cld_all~MCD09_mean|month,data=ts,type="p",auto.key=T)

