## Long term summaries
seasconc <- function(x, na.rm=F) {
  #################################################################################################
  ## Precipitation Concentration function
  ## This function calculates Precipitation Concentration based on Markham's (1970) technique as described in Schulze (1997)
  ## South Africa Atlas of Agrohydology and Climatology - R E Schulze, M Maharaj, S D Lynch, B J Howe, and B Melvile-Thomson
  ## Pages 37-38
  #################################################################################################
  ## x is a vector of monthly quantities
  #################################################################################################
  theta=seq(30,360,30)*(pi/180)                                       # set up angles for each month & convert to radians
  if(!na.rm) if(any(is.na(x))) return(NA)
    rt=sqrt(sum(x * cos(theta),na.rm=T)^2 + sum(x * sin(theta),na.rm=T)^2)    # the magnitude of the summation
    rsum=sum(x,na.rm=T)
    Pc=(1000*rt/rsum)
    Pc
}


## Long term summaries
seastheta <- function(x, na.rm=F) {
  #################################################################################################
  ## Precipitation Concentration function
  ## This function calculates Precipitation Concentration based on Markham's (1970) technique as described in Schulze (1997)
  ## South Africa Atlas of Agrohydology and Climatology - R E Schulze, M Maharaj, S D Lynch, B J Howe, and B Melvile-Thomson
  ## Pages 37-38
  #################################################################################################
  ## x is a vector of 12 monthly quantities
  #################################################################################################
  theta=seq(30,360,30)*(pi/180)                                       # set up angles for each month & convert to radians
  if(!na.rm) if(any(is.na(x))) return(NA)
  xsin=sum(x*sin(theta),na.rm=T)
  xcos=sum(x*cos(theta),na.rm=T)
  xatan=atan(xsin/xcos)
  if(xsin>=0 & xcos>=0)  thetat=abs((180/pi)*xatan)
  if(xsin>=0 & xcos<=0)  thetat=180-abs((180/pi)*xatan)
  if(xsin<=0 & xcos<=0)  thetat=180+abs((180/pi)*xatan)
  if(xsin<=0 & xcos>=0)  thetat=360-abs((180/pi)*xatan)
  thetat=thetat*10
  return(thetat)
}
