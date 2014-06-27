## Assessment of trends in cloud cover
source("analysis/setup.R")

jun=stack("/Users/adamw/GoogleDrive/Work/Cloud/EarthEngineOutput/ee_mcd09cf/20140618_trend_b2599552_20002014_MOD09GA_6-0000000000-0000000000.tif")
dec=stack("/Users/adamw/GoogleDrive/Work/Cloud/EarthEngineOutput/ee_mcd09cf/20140618_trend_b2599552_20002014_MOD09GA_12-0000000000-0000000000.tif")


trends=stack(dec[[1]],jun[[1]])
names(trends)=c("December","June")
NAvalue(trends)=-32767
gain(trends)=.000001  #convert back to change per year
n=100

cellStats(trends,range)
histogram(trends);abline(v=0)

tcoast=spTransform(coast,CRS(projection(trends)))

png("manuscript/figures/trends.png",width=2100,height=2000,res=300,pointsize=42,bg="white")
## set plotting parameters
my.theme = trellis.par.get()
my.theme$strip.background=list(col="transparent")
trellis.par.set(my.theme)

levelplot(trends,col.regions=bgr(values(trends))$col,at=seq(-.030,.030,len=n),maxpixels=1e7)+
  layer(sp.lines(tcoast,col="black",lwd=.5),under=F)

dev.off()
