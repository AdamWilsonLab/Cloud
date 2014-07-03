extract.lm <- function(model) {
  s <- summary(model)
  names <- rownames(s$coef)
  co <- s$coef[, 1]
  se <- s$coef[, 2]
  pval <- s$coef[, 4]
  rs <- s$r.squared
  n <- as.integer(nobs(model))
  rmse=sqrt(mean((residuals(s)^2)))
  gof <- c(rs, rmse, n)
  gof.names <- c("R-Squared","RMSE","n")
  tr <- createTexreg(coef.names = names, coef = co, se = se, 
                     pvalues = pval, gof.names = gof.names, gof = gof)
  return(tr) 
}
setMethod("extract", signature = className("lm", "stats"),definition = extract.lm)
