# Load the libraries and set working directory
source("analysis/setup.R")

library(dismo)
library(grid) # needed for arrow function
library(ggplot2)


#### load summary tables
ds=do.call(rbind.data.frame,lapply(list.files("output/sdm/",pattern="summary.csv",full=T,recursive=T)[2:3],read.csv))
ds$species=sub("_"," ",ds$species)
print(xtable(ds,caption="Evaluation of distribution models using interpolated precipitation or cloud product"), 
      include.rownames = FALSE,"html",file="manuscript/modeltable.html",digits=2)


## Make a plot to explore the data
pdf(file=paste0("manuscript/figures/SDM_",paste(sp,collapse="_"),".pdf"),width=11,height=7)

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
