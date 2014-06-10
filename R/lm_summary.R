lm_summary=function(x,y,var) {
 m=lm(x~y)
 s=summary(m)
 rs <- s$r.squared
 n <- as.integer(nobs(m))
 rmse=sqrt(mean((residuals(s)^2)))
 me=mean(residuals(s),na.rm=T)
 list(n=n,rs=rs,rmse=rmse,me=me)[[var]]
}