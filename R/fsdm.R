## function to fit SDM and write out results and evaluations as a netcdf file to enable better archiving and faster post-processing

fsdm<-function(fdata,data,model,burnin=10000,mcmc=10000,thin=10,nchains=3,...){
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
      burnin=burnin, mcmc=mcmc, thin=thin,
      beta.start=0,
      suitability.pred=data,
      mubeta=0, Vbeta=1.0E6,
      save.p=0,
      verbose=1,
      seed=round(runif(1,0,1000)))
  }

  #hSDM.nc(results=tres1,species=sp2,data=data,fdata=fdata,modelname=modelname,outputdir=outputdir,verbose=T,autocor=T) 
  hSDM.nc(results=tres1,data=data,fdata=fdata,model=model,...) 
  
}


