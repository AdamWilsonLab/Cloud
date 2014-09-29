## function to fit SDM and write out results and evaluations as a netcdf file to enable better archiving and faster post-processing

fsdm<-function(fdata,data,model="~PPTJAN+PPTJUL+PPTSEAS+MAT",modelname="Precipitation",nchains=3,autocor=T,evaluate=T,keepall=F){
 
  require(ncdf4)
  tres1=foreach(ch=1:nchains) %dopar% {
    tres2=hSDM.ZIB(
      suitability=as.character(model),
      presences=fdata$presences,
      observability=~1,
      mugamma=0, Vgamma=1.0E6,
      gamma.start=0,
      trials=fdata$trials,
      data=fdata,
      burnin=1000, mcmc=1000, thin=10,
      beta.start=0,
      suitability.pred=data,
      mubeta=0, Vbeta=1.0E6,
      save.p=0,
      verbose=1,
      seed=round(runif(1,0,1000)))
  }

  ## assess convergence
  coef_list=as.mcmc.list(lapply(tres1,FUN=function(x) x$mcmc))
  c1=gelman.diag(coef_list, confidence = 0.95, transform=FALSE, autoburnin=TRUE, multivariate=TRUE)
  colnames(c1$psrf)=c("GelmanPSRF","GelmanPSRF.CI")
  c2=geweke.diag(coef_list[[1]], frac1=0.1, frac2=0.5)
  c3=heidel.diag(coef_list[[1]], eps=0.1, pvalue=0.05)
  colnames(c3)=c("Heidel.Stationarity","Heidel.Start","Heidel.P","Heidel.Halfwidth","Heidel.Mean","Heidel.Halfwidth")
  
  conv=rbind("Global"=c(c1$mpsrf,rep(NA,8)),
    cbind(c1$psrf,GewekeZ=c2$z,c3))
  
  ## combine the chains for each model
  coef=data.frame(param=colnames(tres1[[1]][[1]]),
                  summary(as.mcmc.list(lapply(tres1,FUN=function(x) x$mcmc)))$statistics,
                  summary(as.mcmc.list(lapply(tres1,FUN=function(x) x$mcmc)))$quantiles)
  pred=data.frame(x=data$x,y=data$y,
                  pred=rowMeans(do.call(cbind,lapply(tres1,FUN=function(x) x$prob.p.pred))),
                  cell=cellFromXY(senv,xy=pred[,c("x","y")]))

  
  predr=rasterFromXYZ(xyz=pred[,c("x","y","pred")])
  projection(predr)='+proj=longlat'
  #writeRaster(predr*100,file=paste0(outputdir,sp2,".tif"),overwrite=T)


##### Autocorrelation
ac1=acor_table(predr[[1]],verbose=T)
ac1$model=names(predr)[1]
ac2=acor_table(predr[[2]],verbose=T)
ac2$model=names(predr)[2]

ac=rbind.data.frame(ac1,ac2)


## Global Spatial Autocorrelation
fspac=paste0(outputdir,sp2,".csv")
if(!file.exists(fspac)){
  spac=do.call(rbind.data.frame,
               lapply(1:nlayers(predr),function(i) {
                 data.frame(species=paste(sp,collapse=" "),
                            model=names(predr)[i],
                            MoransI=Moran(predr[[i]],w=matrix(1,51,51)),
                            GearyC=Geary(predr[[i]],w=matrix(1,51,51)))}
               ))
  write.csv(spac,file=fspac,row.names=F)
}
spac=read.csv(fspac)

## Range size
rs=cellStats(predr,"sum")
(rs[[1]]-rs[[2]])/rs[[1]]

# http://www.bayesian-inference.com/modelfit
# http://voteview.com/DIC-slides.pdf
dic=coef[grepl("Deviance",coef$param),]
dic$pv=(dic$SD^2)/2
dic$dic=dic$Mean+dic$pv

## AUC
aucdat=merge(fdata[,c("presences","trials","cell")],pred,by=c("cell"))

e=evaluate(p=aucdat$pred[aucdat$presences>0],a=aucdat$pred[aucdat$presences==0])

restable=data.frame(nPresence=e@np,nTrials=e@na,auc=e@auc,cor=e@cor)
restable$species=sp2
restable$model=dic$model
restable$Deviance=dic$Mean
restable$Pv=dic$pv
restable$DIC=dic$dic
##
restable[,c("MoransI","GearyC")]=spac[match(restable$model,spac$model),c("MoransI","GearyC")]

restable=restable[order(restable$model),
                  c("species","model","nPresence","nTrials","Deviance","Pv","DIC","auc","cor","MoransI","GearyC")]

############################################################
## Write results to file
  ## Set dimentions
  d_lat=ncdim_def("lat",units="degrees_north",vals=sort(unique(pred$y),decreasing=F),longname="latitude")
  d_lon=ncdim_def("lon",units="degrees_east",vals=sort(unique(pred$x)),longname="longitude")
  d_iter=ncdim_def("iter",units="iterations",vals=1:nrow(tres1[[1]]$mcmc),unlim=FALSE)
  d_params=ncdim_def("parameters",units="",create_dimvar=F,vals=1:nrow(coef),unlim=FALSE)  #posterior parameter summaries
  d_params2=ncdim_def("parameters2",units="",create_dimvar=F,vals=1:ncol(coef),unlim=FALSE)  #posterior parameter summaries
  
  d_converge1=ncdim_def("convergence",units="",create_dimvar=F,vals=1:ncol(conv),unlim=FALSE)  #convergence summaries
  d_converge2=ncdim_def("convergence2",units="",create_dimvar=F,vals=1:nrow(conv),unlim=FALSE)  #convergence summaries
  
  d_eval1=ncdim_def("evaluation",units="",create_dimvar=F,vals=1:ncol(restable),unlim=FALSE)  #evaluation summaries
  d_eval2=ncdim_def("evaluation2",units="",create_dimvar=F,vals=1:nrow(restable),unlim=FALSE)  #evaluation summaries
  
  ## Define variables
  comp=9  #define compression level: 9 is the highest
  v_var_mean=ncvar_def("p",units="Probability",dim=list(d_lon,d_lat),missval=-999,
                       longname="Probability of Occurrence",compress=comp)
  v_var_parameters=ncvar_def("parameters",units="values",dim=list(d_params,d_params2),missval=-999,
                           longname="Convergence Metrics",compress=comp)
    v_var_converge=ncvar_def("convergence",units="values",dim=list(d_converge1,d_converge2),missval=-999,
                        longname="Convergence Metrics",compress=comp)
  v_var_evaluate=ncvar_def("evaluation",units="values",dim=list(d_eval1,d_eval2),missval=-999,
                    longname="Evaluation Metrics",compress=comp)
  
  if(keepall){
    v_var=ncvar_def("p_sample",units="Probability",dim=list(d_lon,d_lat),missval=-999,
                    longname="Probability of Occurrence",compress=comp)
  }
  
  ## set up nc file
  ncfile=paste0(outputdir,"/",sp2,"_",gsub("~|[+]|_","_",model),".nc",sep="")
  if(file.exists(ncfile)) file.remove(ncfile)
  if(keepall) nc_create(ncfile,vars=list(v_var,v_var_mean,v_var_sd,v_krige,v_var_valid),verbose=F)   #save every iteration
  if(!keepall) nc_create(ncfile,vars=list(v_var_mean,v_var_parameters,v_var_converge,v_var_evaluate),verbose=F) 
  
  nc=nc_open(ncfile,write=T)
  print("NetCDF file created, adding data")
  
  ncatt_put(nc,"parameters", "names", paste(colnames(tres1$mcmc[[1]]),collapse=","))
  ncatt_put(nc,"evaluation","columns",paste(colnames(attr(kc,"meta")$valid),collapse=","),prec="char")
  
  ## Add data
  ncvar_put(nc,paste(vari,"_mean",sep=""),vals=var_mean[,ncol(var_mean):1],start=c(1,1,1),c(-1,-1,1),verb=F)
  ncvar_put(nc,paste(vari,"_sd",sep=""),vals=var_sd[,ncol(var_sd):1],start=c(1,1,1),c(-1,-1,1),verb=F)
  ncvar_put(nc,paste(vari,"_valid",sep=""),vals=t(as.matrix(attr(kc,"meta")$valid)),start=c(1,1,1),c(-1,-1,1),verb=F)
  #ncvar_put(nc,paste(vari,"_params",sep=""),vals=apply(kc$posterior$sample,2,function(x) c(mean=mean(x),var=var(x))),start=c(1,1,1),c(-1,-1,1),verb=F)
  ncvar_put(nc,paste(vari,"_params",sep=""),vals=t(as.matrix(kc$posterior$sample)),start=c(1,1,1),c(-1,-1,1),verb=F)

  if(keepall){
    ncvar_put(nc,vari,vals=var_iter,start=c(1,1,1,1),c(-1,-1,-1,1),verb=F)
    ncvar_put(nc,paste(vari,"_q2.5",sep=""),vals=var_q2.5[,ncol(var_q2.5):1],start=c(1,1,1),c(-1,-1,1),verb=F)
    ncvar_put(nc,paste(vari,"_q97.5",sep=""),vals=var_q97.5[,ncol(var_q97.5):1],start=c(1,1,1),c(-1,-1,1),verb=F)
  }
  
  
  print("Data added, updating attributes")
  ################################
  ## Attributes
  ncatt_put(nc,"iter", "description","Posterior samples from predictive distribution",prec="character")
  ## Global Attributes
  ncatt_put(nc,varid=0, "title",paste0("Predicted p(occurrence) for ",sp2),prec="character")
  ncatt_put(nc,varid=0, "institution","Map of Life, Yale University, New Haven, CT",prec="character")
  ncatt_put(nc,varid=0, "source","Modeled Species Distributions",prec="character")
  ncatt_put(nc,varid=0, "comment","Adam M. Wilson (adam.wilson@yale.edu)",prec="character")
  ncatt_put(nc,varid=0, "RunningTime",paste(round(attr(kc,"meta")$duration,2),attr(attr(kc,"meta")$duration,"units")),prec="character")
  
  ## Close the file
  nc_sync(nc)
  nc_close(nc)
  
  
  
  
  
  
}