# Load the libraries and set working directory
source("analysis/setup.R")

library(dismo)

#################
load(fres)

coef=do.call(rbind,lapply(res,function(x) x$coef))
pred=do.call(rbind,lapply(res,function(x) x$pred))
rho=do.call(rbind,lapply(res,function(x) x$rho))

pred=pred[pred$model%in%mods$model[1:2],]
pred$cell=cellFromXY(senv,xy=pred[,c("x","y")])

# http://www.bayesian-inference.com/modelfit
# http://voteview.com/DIC-slides.pdf
dic=coef[grepl("Deviance",coef$param),]
dic$pv=(dic$SD^2)/2
dic$dic=dic$Mean+dic$pv

dic

## AUC
aucdat=merge(fdata[,c("presences","trials","cell")],pred,by=c("cell"))

feval=function(x){
  e=evaluate(p=x$pred[x$presences>0],a=x$pred[x$presences==0])
  data.frame(nPresence=e@np,nTrials=e@na,auc=e@auc,cor=e@cor)
}

restable=do.call(rbind.data.frame,by(aucdat,INDICES=aucdat$modelname,FUN=feval))
restable$species=paste(sp,collapse=" ")
restable$model=unique(as.character(aucdat$modelname))
restable$Deviance=dic$Mean[match(restable$model,dic$modelname)]
restable$Pv=dic$pv[match(restable$model,dic$modelname)]
restable$DIC=dic$dic[match(restable$model,dic$modelname)]

restable=restable[,c("species","model","nPresence","nTrials","Deviance","Pv","DIC","auc","cor")]

write.csv(restable,file=paste0("data/SDM/",paste(sp,collapse="_"),"_summary.csv"),row.names=F)


#### load summaries
ds=do.call(rbind.data.frame,lapply(list.files("data/SDM/",pattern="summary.csv",full=T),read.csv))

print(xtable(ds,caption="Evaluation of distribution models using interpolated precipitation or cloud product"), 
      include.rownames = FALSE,"html",file="manuscript/modeltable.html")

library(grid) # needed for arrow function
library(ggplot2)

## Make a plot to explore the data
pdf(file=paste0("manuscript/figures/SDM_",paste(sp,collapse="_"),".pdf"),width=11,height=7)

gplot(senv) + geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradientn(colours=c('white','blue','red'),na.value="transparent") +
  coord_equal()

ggplot(pred) + geom_tile(aes(x=x,y=y,fill = pred)) +
  facet_wrap(~modelname,ncol=1) +
  scale_fill_gradientn(colours=c('white','blue','red')) +
  geom_point(data=spd_sp[spd_sp$presence==0,]@data,aes(x=lon,y=lat),pch="-",col="grey")+
  geom_point(data=spd_sp[spd_sp$presence==1,]@data,aes(x=lon,y=lat),pch="+")+
  #  xlim(c(17.8,25.2))+ylim(-35,-32)+
  coord_equal()

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


## Fit a simple GLM to the data
obs=as.matrix(data[,c("presences","trials")])
m_glm=lapply(1:nrow(mods),function(i) glm(paste0("obs",mods$formula[i]),data=data,family="binomial"))

screenreg(m_glm,custom.model.names=as.character(mods$name),format="markdown")

#= Parameter estimates
summary(res[[1]])
xyplot(mod1$mcmc)
#= Predictions
summary(mod1$theta.latent)
summary(mod1$theta.pred)
plot(theta,mod1$theta.pred)

### Check Convergence
xyplot(mod1$mcmc, main="Beta",strip=strip.custom(factor.levels=c("intercept",vars)))
gelman.diag(mod1$mcmc,confidence = 0.95,autoburnin=F,multivariate=T)


