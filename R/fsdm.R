## function to fit SDM and write out results and evaluations as a netcdf file to enable better archiving and faster post-processing

P.hSDM.ZIB<-function(nchains=3,...){
  foreach(ch=1:nchains) %dopar% {
    hSDM.ZIB(...)
  }
}


