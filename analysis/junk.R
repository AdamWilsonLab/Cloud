### identify elevation band
adata=data.frame(values(stack(env[["ALT"]],spr)))
adata=na.omit(adata)
adata=adata[adata$trials>0,]
adata$p=adata$presences>0

ggplot(adata, aes(x=ALT, fill=p)) + geom_density(alpha=.3)+
  geom_vline(aes(xintercept=c(900,2400),col="red",lwd=2))+
  geom_rug(aes(group=p),sides="b") 

ggplot(adata, aes(y=ALT, x=p)) + 
  geom_boxplot()+ geom_jitter()+
  geom_hline(aes(yintercept=c(900,2400),col="red",lwd=2))+
  geom_rug(aes(group=p),sides="b") 

altrange=calc(env[["ALT"]],function(x) x>900&x<2400)
names(altrange)="altrange"
#senv=stack(senv,altrange)


#splom(senv)



## Correlogram of points
#library(ncf)
#s=sample(1:nrow(fdata),100)
#pcg=correlog(x=fdata$x[s],y=fdata$y[s],z=fdata$presence[s]>0,latlon=T,increment=50)

#library(spdep)
#ds=c(seq(0,100,by=5))
#joincount.mc(HICRIME, nb2listw(COL.nb, style="B"), nsim=99)
#d1=dnearneigh(as.matrix(fdata[,c("x","y")]),d1=0,d2=100,longlat=T)
#plot(pcg);abline(h=0,col="red")

#colnames(coef)[grepl("X",colnames(coef))]=paste0("Q",sub("X","",colnames(coef)[grepl("X",colnames(coef))]))
#res$param=sub("beta[.]","",rownames(res))
#0resl=melt(res,id.vars=c("model","suitability","param"))
#levels(resl$variable)=c(0.025,0.25,0.5,0.75,0.975)
library(spdep)
set.ZeroPolicyOption(T)  

nb=dnearneigh(spd,0,5,longlat=T)
plot(nb)
nblist=nb2listw(nb,style="B")
#tdrop=do.call(c,lapply(nblist$weights,function(x) length(x)==0))
#nblist=nb2listw(nb,style="B",zero.policy=TRUE)
joincount.test(as.factor(spd$presence>0),nblist)






penv=gplot(senv) + geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colours=c('white','blue','red'),na.value="transparent") +
  coord_equal()+ theme(legend.position="bottom")

## compare predictions
p1=
  gplot(predr,maxpixels=1e4) + geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colours=c('white','blue','red'),values = c(0,.6,.75,1), 
                       rescaler = function(x, ...) x,
                       na.value="transparent") +
  coord_equal()+xlab("Longitude")+ylab("Latitude")+
  geom_point(data=spd[spd$presence==0,]@data,aes(x=lon,y=lat),pch=1,col=grey(.8),cex=1)+
  geom_point(data=spd[spd$presence==1,]@data,aes(x=lon,y=lat),pch=3,cex=2)+
  #  geom_line(data=fcoast,mapping=aes(long,lat,group=group)) +
  labs(fill = "p(presence)")+ 
  theme(legend.position="bottom")

p2=ggplot(aucdat, aes(factor(presences>0,labels=c("Absent","Present")), pred))+ 
  geom_boxplot()+facet_wrap(~model)+
  xlab("Observed")+ylab("p(presence)")+
  coord_flip()

## plot the correlogram
p3=ggplot(ac, aes(x=dist2, y=mean, group=model))+
  geom_line(aes(colour=model))+
  geom_pointrange(aes(ymax = max, ymin=min, group=model,colour=model))+
  geom_line(col="black")+
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  ylab("Autocorrelation")+xlab("Distance (km)")+
  theme(legend.position=c(.25, .25))

p4=  gplot(predr,maxpixels=1e7) + geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colours=c('white','blue','red'),values = c(0,.6,.75,1), 
                       rescaler = function(x, ...) x,
                       na.value="white") +
  coord_equal()+xlab("Longitude")+ylab("Latitude")+
  geom_point(data=spd[spd$presence==0,]@data,aes(x=lon,y=lat),pch=1,col=grey(.8),cex=1)+
  geom_point(data=spd[spd$presence==1,]@data,aes(x=lon,y=lat),pch=3,cex=2)+
  labs(fill = "p(presence)")+ 
  theme(legend.position="bottom")+
  ylim(c(-1,7.5))+xlim(c(-79.5,-72))

## make the pdf
pdf(file=paste0(outputdir,"env_",sp2,".pdf"),width=5,height=10)
penv
dev.off()


pdf(file=paste0(outputdir,"SDM_",sp2,".pdf"),width=5,height=10)
grid.newpage()

vp1 <- viewport(width = 1, height = 0.3, x=.5,y=.15)
vp2 <- viewport(width = 1, height = .5, x = .5, y = 0.55)
vp3 <- viewport(width = .95, height = .2, x = .5, y = 0.9)

print(p3, vp = vp1)
print(p4, vp = vp2)
print(p2, vp = vp3)

dev.off()



ggplot(coef[coef$param!="Deviance",], aes(x = Mean, y = param,colour=model))+  
  geom_segment(aes(xend = X97.5.,yend=param),lwd=2)


## regional
#p1=ggplot(pred) + geom_tile(aes(x=x,y=y,fill = pred)) +
#  facet_wrap(~modelname,ncol=1) +
#  scale_fill_gradientn(colours=c('white','blue','red'),
#                       name="P(Presence)") +
#  coord_equal()+
#  xlim(c(20,25))+ylim(-34.2,-33.5)+
#  geom_point(data=prot[prot$pro!=sp,]@data,aes(x=londd,y=latdd),col="grey",alpha=.2,pch=1)+
#  geom_point(data=prot[prot$pro==sp,]@data,aes(x=londd,y=latdd),pch="+",cex=1)+
#  theme(panel.background = element_rect(fill='transparent'),legend.key.width = unit(1.5, "cm"))+
#  xlab(label="Longitude")+ylab("Latitude")

# p2=gplot(env[["alt"]],maxpixels=5e6) + geom_tile(aes(fill = value)) +
#   facet_grid(~variable,labeller=function(...) return("Protea cynaroides occurrence")) +
#     scale_fill_gradientn(colours=c('darkgreen','yellow','red'),na.value="transparent",
#                        name="Elevation (m)") +
#   coord_equal()+
#   geom_point(data=prot[prot$pro!=sp,]@data,aes(x=londd,y=latdd),col="grey",alpha=.2,pch=1,cex=.5)+
#   geom_point(data=prot[prot$pro==sp,]@data,aes(x=londd,y=latdd),pch="+",cex=2)+
#   scale_colour_manual(name = 'Occurrence', 
#                       values =c('black'='black','grey'='grey'),
#                       labels = c('Presence','Absence'))+
#   xlim(c(20,25))+ylim(-34.2,-33.5)+
#   xlab(label="")+
#   theme(panel.background = element_rect(fill='transparent'),legend.key.width = unit(1.5, "cm"))+
#   ylab("")


# grid.newpage()
# vp1 <- viewport(width = 1, height = 0.67, x=.5,y=.33)
# vp2 <- viewport(width = 1, height = .33, x = .5, y = 0.72)
# print(p1, vp = vp1)
# print(p2, vp = vp2)
dev.off()

