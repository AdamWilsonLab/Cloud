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



extract.glm <- function (model, include.sdm=TRUE, include.aic = TRUE, include.bic = TRUE, include.loglik = TRUE, 
            include.deviance = TRUE, include.nobs = TRUE, ...) 
  {
    s <- summary(model, ...)
    coefficient.names <- rownames(s$coef)
    coefficients <- s$coef[, 1]
    standard.errors <- s$coef[, 2]
    significance <- s$coef[, 4]
    aic <- AIC(model)
    bic <- BIC(model)
    lik <- logLik(model)[1]
    dev <- deviance(model)
    n <- nobs(model)
    gof <- numeric()
    gof.names <- character()
    gof.decimal <- logical()
    ## add some SDM things
    if (include.sdm == TRUE) {
    response=as.character(model$terms)[2]
    e1=evaluate(model=model,p=model$fitted.values[model$model[,response]==1],a=model$fitted.values[model$model[,response]==0])
    auc=e1@auc
    cor=e1@cor
    equal_sens_spec=threshold(e1)$equal_sens_spec
    gof <- c(gof, auc,cor)#,equal_sens_spec)
    gof.names <- c(gof.names, "AUC","COR")#,"Threshold (sensitivity=specificity)")
    gof.decimal <- c(gof.decimal, TRUE, T)#, T)
    }
    ###   
    if (include.aic == TRUE) {
      gof <- c(gof, aic)
      gof.names <- c(gof.names, "AIC")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    if (include.bic == TRUE) {
      gof <- c(gof, bic)
      gof.names <- c(gof.names, "BIC")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    if (include.loglik == TRUE) {
      gof <- c(gof, lik)
      gof.names <- c(gof.names, "Log Likelihood")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    if (include.deviance == TRUE) {
      gof <- c(gof, dev)
      gof.names <- c(gof.names, "Deviance")
      gof.decimal <- c(gof.decimal, TRUE)
    }
    if (include.nobs == TRUE) {
      gof <- c(gof, n)
      gof.names <- c(gof.names, "Num. obs.")
      gof.decimal <- c(gof.decimal, FALSE)
    }
    tr <- createTexreg(coef.names = coefficient.names, coef = coefficients, 
                       se = standard.errors, pvalues = significance, gof.names = gof.names, 
                       gof = gof, gof.decimal = gof.decimal)
    return(tr)
  }

  
setMethod("extract", signature = className("glm", "stats"),definition = extract.glm)
