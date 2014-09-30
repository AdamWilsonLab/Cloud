hSDM.nc<-function(results,species,modelname,data,fdata,outputdir,autocor=T,keepall=F,verbose=T){

  today=format(Sys.Date(),format="%Y%m%d")
  time=format(Sys.time(),format="%H%M%S")
  ncfile=paste0(outputdir,"/",species,"_",modelname,"_",today,"_",time,".nc",sep="")
  if(verbose) writeLines(paste("Preparing to write",ncfile))
  
  ## assess convergence for each parameter
  if(verbose) writeLines("Calculating convergence metrics")
  parameters_list=mcmc.list(lapply(results,FUN=function(x) x$mcmc))
  c1=gelman.diag(parameters_list, confidence = 0.95, transform=FALSE, autoburnin=FALSE, multivariate=TRUE)
  colnames(c1$psrf)=c("GelmanPSRF","GelmanPSRF.CI")
  c2=geweke.diag(parameters_list[[1]], frac1=0.1, frac2=0.5)
  c3=heidel.diag(parameters_list[[1]], eps=0.1, pvalue=0.05)
  colnames(c3)=c("Heidel.Stationarity","Heidel.Start","Heidel.P","Heidel.HalfwidthTest","Heidel.Mean","Heidel.Halfwidth")
  class(c3)="matrix"  #convert class for easier cbinding
  
  ## summarize posterior distributions  
  if(verbose) writeLines("Calculating posterior summaries")
  parameters=data.frame(mean=summary(parameters_list)$statistics[,"Mean"],
                        sd=summary(parameters_list)$statistics[,"SD"],
                        median=summary(parameters_list)$quantiles[,"50%"],
                        HPDinterval(mcmc(as.matrix(parameters_list))),
                        RejectionRate=rejectionRate(parameters_list),
                        c1$psrf,GewekeZ=c2$z,c3)
  
  ## predictions for each cell
  if(verbose) writeLines("Summarizing pixel-level posteriors")
  
  pred=data.frame(x=data$x,y=data$y,
                  pred=rowMeans(do.call(cbind,lapply(results,FUN=function(x) x$prob.p.pred))),
                  cell=cellFromXY(senv,xy=pred[,c("x","y")]))
  predr=rasterFromXYZ(xyz=pred[,c("x","y","pred")])
  projection(predr)='+proj=longlat'
  
  ##### Autocorrelation
  if(autocor){
    if(verbose) writeLines("Calculating autocorrelation of output")
    ac=acor_table(predr,verbose=F)
    ## Global Spatial Autocorrelation
    spac=c(MoransI=Moran(predr,w=matrix(1,11,11)),
           GearyC=Geary(predr,w=matrix(1,11,11)))
  }
  
  
  ## AUC
  if(verbose) writeLines("Calculating model evaluation metrics")
    aucdat=merge(fdata[,c("presences","trials","cell")],pred,by=c("cell"))
  e=evaluate(p=aucdat$pred[aucdat$presences>0],a=aucdat$pred[aucdat$presences==0])
  
  evaluation=data.frame(nPresence=e@np,nTrials=e@na,auc=e@auc,cor=e@cor)
  evaluation$Deviance=parameters$mean[grepl("Deviance",rownames(parameters))]
  evaluation$Pd=(parameters$sd[grepl("Deviance",rownames(parameters))]^2)/2 #Gelman BDA pg 182
  evaluation$DIC=evaluation$Deviance+evaluation$Pd
  evaluation$nChains=nchains
  evaluation$nBurnin=mcpar(tres1[[1]]$mcmc)[1]
  evaluation$nIter=mcpar(tres1[[1]]$mcmc)[2]
  evaluation$thin=mcpar(tres1[[1]]$mcmc)[3]
  evaluation$GlobalGelmanMPSRF=c1$mpsrf
  if(autocor) evaluation[,c("MoransI","GearyC")]=spac
  
  ############################################################
  ## Write results to netcdf file
  if(verbose) writeLines("Setting up netCDF file")
  
  ## Set dimentions
  d_lat=ncdim_def("lat",units="degrees_north",vals=sort(unique(pred$y),decreasing=F),longname="latitude")
  d_lon=ncdim_def("lon",units="degrees_east",vals=sort(unique(pred$x)),longname="longitude")
  d_iter=ncdim_def("iter",units="iterations",longname="Posterior Iterations",vals=1:nrow(results[[1]]$mcmc),unlim=TRUE)
  
  d_params1=ncdim_def("parameters1",units="",create_dimvar=F,vals=1:nrow(parameters),unlim=FALSE)  #posterior parameter summaries
  d_params2=ncdim_def("parameters2",units="",create_dimvar=F,vals=1:ncol(parameters),unlim=FALSE)  #posterior parameter summaries
  
  d_eval1=ncdim_def("evaluation1",units="",create_dimvar=F,vals=1:ncol(evaluation),unlim=FALSE)  #evaluation summaries
  d_eval2=ncdim_def("evaluation2",units="",create_dimvar=F,vals=1:nrow(evaluation),unlim=FALSE)  #evaluation summaries
  
  if(autocor){
    d_ac1=ncdim_def("autocorrelation1",units="",create_dimvar=F,vals=1:ncol(ac),unlim=FALSE)  #autocorrelation summaries
    d_ac2=ncdim_def("autocorrelation2",units="",create_dimvar=F,vals=1:nrow(ac),unlim=FALSE)  #autocorrelation summaries
  }
  
  ## Define variables
  comp=9  #define compression level: 9 is the highest
  v_var_mean=ncvar_def("p",units="Probability",dim=list(d_lon,d_lat),missval=-999,
                       longname="Probability of Occurrence *100",compress=comp,prec="integer")
  v_var_parameters=ncvar_def("parameters",units="values",dim=list(d_params1,d_params2),missval=-999,
                             longname="Posterior Parameter Values",compress=comp)
  v_var_evaluate=ncvar_def("evaluation",units="values",dim=list(d_eval1,d_eval2),missval=-999,
                           longname="Evaluation Metrics",compress=comp)
  
  if(autocor){
    v_var_ac=ncvar_def("ac",units="values",dim=list(d_ac1,d_ac2),missval=-999,
                       longname="Autocorrelation",compress=comp)
  }
  ###  if(keepall){
  #    v_var=ncvar_def("p_sample",units="Probability",dim=list(d_lon,d_lat),missval=-999,
  #                    longname="Probability of Occurrence",compress=comp)
  #  }
  
  ## set up nc file
  if(file.exists(ncfile)) file.remove(ncfile)
  if(keepall) nc_create(ncfile,vars=list(v_var,v_var_mean,v_var_sd,v_krige,v_var_valid),verbose=F)   #save every iteration
  if(!keepall) nc_create(ncfile,vars=list(v_var_mean,v_var_parameters,v_var_evaluate,v_var_ac),verbose=F) 
  
  nc=nc_open(ncfile,write=T)
  if(verbose) print("NetCDF file created, adding data")
  
  ncatt_put(nc,"parameters", "colnames", paste(parameters$param,collapse=","))
  ncatt_put(nc,"parameters", "rownames", paste(colnames(parameters),collapse=","))
  
  ncatt_put(nc,"evaluation","colnames",paste(colnames(evaluation),collapse=","),prec="char")
  ncatt_put(nc,"evaluation","rownames",paste(rownames(evaluation),collapse=","),prec="char")
  
  if(autocor) ncatt_put(nc,"ac","colnames",paste(colnames(ac),collapse=","),prec="char")
  
  ## Add data
  ncvar_put(nc,"parameters",vals=as.matrix(parameters),start=c(1,1),c(-1,-1),verb=F)
  ncvar_put(nc,"evaluation",vals=as.matrix(evaluation),start=c(1,1),c(-1,-1),verb=F)
  if(autocor) ncvar_put(nc,"ac",vals=as.matrix(ac),start=c(1,1),c(-1,-1),verb=F)
  
  ## Add map data
  ncvar_put(nc,"p",vals=1000*t(as.matrix(predr))[,nrow(predr):1],start=c(1,1),c(-1,-1),verb=F)
  ncatt_put(nc,varid="p", "projection",projection(predr),prec="character")
  ncatt_put(nc,varid="p", "projection_format","PROJ.4",prec="character")
  
  
  if(verbose) print("Data added, updating attributes")
  ################################
  ## Attributes
  ## Global Attributes
  ncatt_put(nc,varid=0, "Conventions","Cf-1.4",prec="character")
  ncatt_put(nc,varid=0, "title",paste0("Predicted p(occurrence) for ",sp2),prec="character")
  ncatt_put(nc,varid=0, "institution","Map of Life, Yale University, New Haven, CT",prec="character")
  ncatt_put(nc,varid=0, "source","Modeled Species Distributions",prec="character")
  ncatt_put(nc,varid=0, "comment","Adam M. Wilson (adam.wilson@yale.edu)",prec="character")
  ncatt_put(nc,varid=0, "model",model,prec="character")
  ncatt_put(nc,varid=0, "modelname",modelname,prec="character")
  ncatt_put(nc,varid=0, "species",sp2,prec="character")
  ncatt_put(nc,varid=0, "date",today,prec="character")
  
  ## Close the file
  nc_sync(nc)
  nc_close(nc)
  
  if(verbose) print("Finished...")
  
  #ncvar_put(nc,paste(vari,"_params",sep=""),vals=apply(kc$posterior$sample,2,function(x) c(mean=mean(x),var=var(x))),start=c(1,1,1),c(-1,-1,1),verb=F)
  #ncvar_put(nc,paste(vari,"_params",sep=""),vals=t(as.matrix(kc$posterior$sample)),start=c(1,1,1),c(-1,-1,1),verb=F)
  
  #if(keepall){
  #  ncvar_put(nc,vari,vals=var_iter,start=c(1,1,1,1),c(-1,-1,-1,1),verb=F)
  #  ncvar_put(nc,paste(vari,"_q2.5",sep=""),vals=var_q2.5[,ncol(var_q2.5):1],start=c(1,1,1),c(-1,-1,1),verb=F)
  #  ncvar_put(nc,paste(vari,"_q97.5",sep=""),vals=var_q97.5[,ncol(var_q97.5):1],start=c(1,1,1),c(-1,-1,1),verb=F)
  #}
}