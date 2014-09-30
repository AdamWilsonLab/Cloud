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
