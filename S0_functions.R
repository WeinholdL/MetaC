## this file contains:
#### helper_inverseVar: helper-function to calculate inverse variances
#### returnNA: helper-functio; returns NA in right format
#### withCallingHandlers_LW: customized tryCatch-function for all meta-analysis
####                         and -regression models
#### metamean_raw: meta-analysis (ignoring time-dependency)
#### metareg_raw_tau: meta-regression with covariable tau
#### metamean_logit: meta-analysis on logit-transformed C-index (ignoring time-dependency)
#### metamean_logit_tau: meta-regression on logit-transformed C-index with covariable tau
#### metamean_log: meta-analysis on log-transformed C-index (ignoring time-dependency)
#### metamean_log_tau: meta-regression on log-transformed C-index (ignoring time-dependency)
#### metamean_asin: meta-analysis on asin-transformed C-index (ignoring time-dependency)
#### metamean_asin_tau: meta-regression on asin-transformed C-index (ignoring time-dependency)
#### metareg_tau_spline: meta-regression with restricted cubic spline
#### metareg_tau_mfp: meta-regression with fractional polynomials
#### metamean_MichaelisMenten_nlme: Michealis-Menten meta-regression
#### metamean_expdecay_lw3_R0free: exponential decay meta-regression
#### 
#### 

## libraries:

library(meta)
library(metafor)
library(nlme)
library(rms)
library(pracma)
library(evd)
require(MASS)
library(VGAM)


####----------------------------------------------------------------------------
#### simulate Weibull distributed survival times
getWeibullData <- function(seed, n, gumbelsc=1, gumbelLoc=0, mu=0, sd1=1, sc = 0.5, 
                           beta=2, censrate=2,  silent=T) {
  set.seed(seed)
  
  noise   <- rgumbel(n, loc = gumbelLoc, scale = gumbelsc) # noise term (Gumbel-Distributed)
  beta   <- beta # example: age slope coefficients
  
  X1        <-  matrix(rnorm(n, mean=mu, sd = sd1)) # predictors, e.g. age
  dat        <- data.frame(model.matrix(~ X1)) # data frame
  
  dat$logT  <- as.numeric(- X1*beta + noise*sc) # survival times on log scale -> 
  # Gumbel distributed (Gloc= -X1*beta+noise*sc, Gscale=gumbelsc*sc)
  #   survreg's scale  =    1/(rweibull shape)
  #   survreg's intercept = log(rweibull scale)
  
  R2      <- var(X1%*%beta)/var(dat$logT)
  dat$Ttrue <- exp(- dat$logT) # uncensored survival times, 
  # Weibull distributed with parameter (1/Gscale, exp(-Gloc))
  
  # random censoring times, exponentially distributed,
  # rate = log(2) / median
  # median = log(2)/rate -> rate = log(2)/0.6 -> median = log(2)/rate = 0.6
  dat$Ctrue <- rexp(n, rate = censrate)#log(2)/0.5)#median(dat$Ttrue))
  dat$status <- (dat$Ttrue <= dat$Ctrue) + 0 # ca. 50% censoring
  dat$Ttilde <- pmin(dat$Ttrue, dat$Ctrue)
  dat$id     <- 1:nrow(dat)
  
  return(list("data"=dat, "X"=X1, "beta"=beta, "R2"=R2))
}

####----------------------------------------------------------------------------
#### inverse variances
helper_inverseVar<-function(meta1){
  meta1<-meta1[order(meta1$tau),]
  ## inverse variance
  meta1$invVar<-1/(meta1$SEcindex)^2
  meta1$viweight<-meta1$invVar/sum(meta1$invVar)
  
  ## inverse variance logit-trafo
  meta1$invVar_logit<-1/(meta1$SElogitcindex)^2
  meta1$viweight_logit<-meta1$invVar_logit/sum(meta1$invVar_logit, na.rm=T)
  
  ## inverse variance asin-trafo
  if(!is.null(meta1$SEasincindex)){
    meta1$invVar_asin<-1/(meta1$SEasincindex)^2
    meta1$viweight_asin<-meta1$invVar_asin/sum(meta1$invVar_asin, na.rm=T)
  }
  
  return(meta1)
}



####----------------------------------------------------------------------------
#### initialize result data-frame:

returnNA<-function(name, desc, meta1){
  res<-data.frame("name"=as.character(name),
                  "description" = as.character(desc),
                  "var_method"=factor(NA,levels = c("DL", "REML")),
                  "nstudies"=nrow(meta1),
                  "pred_80"=-9,
                  "se_80"=-9,
                  "CI_80_lower"=-9,
                  "CI_80_upper"=-9,
                  "inCI"=-9,
                  "AbC"= -9,
                  "RMSE"=-9,
                  "Q"=-9, "tau2"=-9, "I2"=-9, "H"=-9,
                  "anyError"= -9, "ErrorMessage"= "",
                  "anyWarning"= -9,"WarningMessage"= "")
  
  return(res)
}
####----------------------------------------------------------------------------
#### customized try-catch for meta-analysis and -regression models

withCallingHandlers_LW<-function(chstr, data){
  anyError <- anyWarning <- FALSE
  err <-warn <-""
  metares<-NULL
  withCallingHandlers(tryCatch(eval(parse(text=chstr)), 
                               error = function(e) {
                                 err <<-e$message
                                 anyError  <<- TRUE}),
                      warning = function(w){
                        warn <<-w$message
                        anyWarning <<-TRUE
                        invokeRestart("muffleWarning")
                      })
  
  WarningError<-data.frame(anyError, err,  anyWarning, warn, stringsAsFactors = FALSE)
  names(WarningError)<-c("anyError",  "ErrorMessage", "anyWarning","WarningMessage")
  
  return(list(metares=metares,WarningError=WarningError))
}


####----------------------------------------------------------------------------
#### time-independent meta-analysis
metamean_raw<-function(meta1, meta1P, co){
  
  #### load true C-index
  load("CIndex_True_Integral_sc1.1_beta0.5_mu0_sd0.5_gumbelLoc0_gumbelScale0.2.rda")
  trueC<-C[C$tau>=min(meta1P$tau) & C$tau<=max(meta1P$tau),]

  #### initialize result data frames:
  pred<-c("DL_raw"=NA, "REML_raw"=NA)
  res<-rbind(returnNA(name="raw", desc="error", meta1=meta1),
             returnNA(name="raw", desc="error", meta1=meta1))
  
  WarningError<-data.frame("var_method"=c(NA,NA),"anyError"=c(NA, NA), 
                           "ErrorMessage"=c(NA, NA), 
                           "anyWarning"=c(NA, NA), "WarningMessage"=c(NA, NA))
  
  
  meta1$SDcindex<-meta1$SEcindex*sqrt(meta1$n)
  
  methods<-c("DL", "REML") ## loop for DL and REML
  for(k in seq_along(methods)){

    #### run metamean function
    chstr<-paste0("metares<-metamean(n = data$n,  mean=data$cindex, sd=data$SDcindex,
                          data = data, method.tau ='",methods[k],"', hakn = TRUE)")
    tempres<-withCallingHandlers_LW(chstr = chstr, data=meta1)
    #### store possible warnings and errors
    WarningError[k,1]<-methods[k]
    WarningError[k,2:5]<-tempres$WarningError    
    #### meta-analysis model
    metares<-tempres$metares
    
    if(WarningError$anyError[k]){
      res[k,]<-returnNA(name="raw", desc="error", meta1=meta1)
      res[k,c("anyError",  "ErrorMessage", "anyWarning","WarningMessage")]<-  WarningError[k,2:5]
    }else{
      #### meta-analysis estimate
      pred[k]<-metares$TE.random
      #### calculate area between prediction and true C-index
      trueC$pred<-metares$TE.random
      AbC<-trapz(x = trueC$tau,y = abs(trueC$pred-trueC$cindex)) / (max(trueC$tau)-min(trueC$tau))
      #### calculate RMSE
      rmse<-sqrt(mean((meta1P$cindex-metares$TE.random)^2))
      
      #### Confidence intervals for estimate
      cilower<-metares$lower.random
      ciupper<-metares$upper.random
      
      #### Coverage
      trueCindex<-C$cindex[C$tau==co]
      covC<-1*((trueCindex>cilower) & (trueCindex<ciupper))
      
      res[k,]<-data.frame("name"="raw",
                          "description" = "",
                          "var_method"=methods[k],
                          "nstudies"=nrow(meta1),
                          "pred_80"=metares$TE.random,
                          "se_80"=metares$seTE.random,
                          "CI_80_lower"=cilower,
                          "CI_80_upper"=ciupper,
                          "inCI"=covC,
                          "AbC"=AbC,
                          "RMSE"=rmse,
                          "Q"=metares$Q, "tau2"=metares$tau2, 
                          "I2"=metares$I2, "H"=metares$H,
                          "anyError"= WarningError[k,2]*1,
                          "ErrorMessage"= WarningError[k,3],
                          "anyWarning"= WarningError[k,4]*1,
                          "WarningMessage"= WarningError[k,5])
    }
  }
  
  reslist<-list(res=res, pred=pred)
  return(reslist)
}


###############################################
metareg_raw_tau<-function(meta1, co){ 
  meta1<-meta1[order(meta1$tau),]
  meta1<-helper_inverseVar(meta1)
  
  #### load true C-index
  load("CIndex_True_Integral_sc1.1_beta0.5_mu0_sd0.5_gumbelLoc0_gumbelScale0.2.rda")
  Ctrue<-C
  Ctrue<-Ctrue[Ctrue$tau>=min(meta1$tau) & Ctrue$tau<=max(meta1$tau),]
  
  #### initialize result data frames:
  nr<-rep(NA, nrow(meta1))
  ppp<-data.frame("pred_DL_raw_tau"=nr, "se_DL_raw_tau"=nr, 
                  "ci.lb_DL_raw_tau"=nr, "ci.ub_DL_raw_tau"=nr,
                  "pi.lb_DL_raw_tau"=nr, "pi.ub_DL_raw_tau"=nr,
                  "pred_REML_raw_tau"=nr, "se_REML_raw_tau"=nr,
                  "ci.lb_REML_raw_tau"=nr, "ci.ub_REML_raw_tau"=nr, 
                  "pi.lb_REML_raw_tau"=nr, "pi.ub_REML_raw_tau"=nr)
  WarningError<-data.frame("var_method"=c(NA,NA),"anyError"=c(NA, NA), 
                           "ErrorMessage"=c(NA, NA), 
                           "anyWarning"=c(NA, NA), "WarningMessage"=c(NA, NA))
  res<-rbind(returnNA(name="raw", desc="error", meta1=meta1),
             returnNA(name="raw", desc="error", meta1=meta1))
  
  
  methods<-c("DL", "REML")
  p1<-1
  for(k in seq_along(methods)){
    #### run meta-regression model
    data <- try(metagen(TE = cindex,
                        seTE = SEcindex,  studlab = id,
                        data = meta1, sm = "RMD", #SMD
                        comb.fixed = TRUE, comb.random = TRUE,
                        title = "C-Index", method.tau=methods[k]))
    chstr<-paste0("metares<-metareg(x = data, formula =~  tau, method='",methods[k],"', hakn = TRUE)")
    tempres<-withCallingHandlers_LW(chstr = chstr, data=data)
    #### store warnings and errors
    WarningError[k,]<-c(methods[k],tempres$WarningError)
    #### meta-analysis model
    metaresRE<-tempres$metares
    
    
    if(WarningError$anyError[k]){
      res[k,]<-returnNA(name="raw_tau", desc="failed", meta1=meta1)
      res[k,c("anyError",  "ErrorMessage", "anyWarning","WarningMessage")]<-  WarningError[k,2:5]
      reslist<-list(res=res, pred=ppp, WarningError=WarningError)
    }else{
      
      #### prediction at 80perc. of max FU:
      newdat<-c("tau"=co)
      pred_80<-as.data.frame(predict(object = metaresRE,newmods = newdat))
      #### confidence intervals of predictions
      cilower<-pred_80$ci.lb
      ciupper<-pred_80$ci.ub  
      
      #### calculate area between prediction and true C-index
      Ctrue$pred<-predict(object = metaresRE,newmods = as.matrix(Ctrue$tau))$pred
      AbC<-trapz(x = Ctrue$tau,y = abs(Ctrue$pred-Ctrue$cindex)) / (max(Ctrue$tau)-min(Ctrue$tau))
      
      #### RMSE
      rmse<-sqrt(mean((meta1$cindex-predict(metaresRE)$pred)^2))
      
      #### Coverage
      trueCindex<-C$cindex[C$tau==co]
      covC<-1*((trueCindex>cilower) & (trueCindex<ciupper))
      
      res[k,]<- data.frame("name"="raw_tau",
                           "description" = "",
                           "var_method"=methods[k],
                           "nstudies"=nrow(meta1),
                           "pred_80"=pred_80$pred,
                           "se_80"=pred_80$se,
                           "CI_80_lower"=cilower,
                           "CI_80_upper"=ciupper,
                           "inCI"=covC,
                           "AbC"=AbC,
                           "RMSE"=rmse,
                           "Q"=metaresRE$QE, 
                           "tau2"=metaresRE$tau2, 
                           "I2"=metaresRE$I2/100, 
                           "H"=metaresRE$H,
                           "anyError"= WarningError[k,2]*1, "ErrorMessage"= WarningError[k,3],
                           "anyWarning"= WarningError[k,4]*1,"WarningMessage"= WarningError[k,5])
      
      ppp1<-as.data.frame(predict(metaresRE))[,1:6]
      names(ppp1)<-paste0(names(ppp1), "_", methods[k])
      ppp[,p1:(p1+5)]<-ppp1
      reslist<-list(res=res, pred=ppp, WarningError=WarningError)
    }
    p1<-p1+6
  }
  
  return(reslist)
  
}



###############################################
metamean_logit<-function(meta1, meta1P, co){  
  
  #### load true C-index
  load("CIndex_True_Integral_sc1.1_beta0.5_mu0_sd0.5_gumbelLoc0_gumbelScale0.2.rda", verbose=T)
  Ctrue<-C
  Ctrue<-Ctrue[Ctrue$tau>=min(meta1P$tau) & Ctrue$tau<=max(meta1P$tau),]
  
  meta1<-meta1[!is.na(meta1$SElogitcindex),]
  meta1<-data.frame("n"=meta1$n, "cindex"=meta1$logitcindex, "sdcindex" = meta1$SElogitcindex*sqrt(meta1$n))
  
  #### initialize result data frames:
  ppp<-c("DL_logit"=NA, "REML_logit"=NA)
  WarningError<-data.frame("var_method"=c(NA,NA),"anyError"=c(NA, NA), "ErrorMessage"=c(NA, NA), 
                           "anyWarning"=c(NA, NA), "WarningMessage"=c(NA, NA))
  res<-rbind(returnNA(name="raw", desc="error", meta1=meta1),
             returnNA(name="raw", desc="error", meta1=meta1))
  
  methods<-c("DL", "REML")
  for(k in seq_along(methods)){
    #### run meta-regression model
    chstr<-paste0("metares<-metamean(n = n,  mean=cindex, sd=sdcindex, data=data, method.tau='",methods[k],"', hakn = TRUE)")
    tempres<-withCallingHandlers_LW(chstr = chstr,data=meta1)
    #### meta-analysis model
    metares<-tempres$metares
    #### store warnings and errors
    WarningError[k,]<-c(methods[k],tempres$WarningError)
    
    if(WarningError$anyError[k]){
      res[k,]<-returnNA(name="logit", desc="failed", meta1=meta1)
      res[k,c("anyError",  "ErrorMessage", "anyWarning","WarningMessage")]<-  WarningError[k,2:5]
      reslist<-list(res=res, pred=ppp, WarningError=WarningError)
    }else{
      
      #### estimate and confidence intervals
      Ctrue$pred<-plogis(metares$TE.random)
      cilower<-plogis(metares$lower.random)
      ciupper<-plogis(metares$upper.random)
      
      #### calculate area between prediction and true C-index
      AbC<-trapz(x = Ctrue$tau,y = abs(Ctrue$pred-Ctrue$cindex)) / (max(Ctrue$tau)-min(Ctrue$tau))
      rmse<-sqrt(mean((meta1P$cindex-plogis(metares$TE.random))^2))
      #### Coverage
      trueCindex<-C$cindex[C$tau==co]
      covC<-1*((trueCindex>cilower) & (trueCindex<ciupper))
      
      res[k,]<-data.frame("name"="logit",
                          "description" = "",
                          "var_method"=methods[k],
                          "nstudies"=nrow(meta1),
                          "pred_80"=plogis(metares$TE.random),
                          "se_80"=metares$seTE.random,
                          "CI_80_lower"=cilower,
                          "CI_80_upper"=ciupper,
                          "inCI"=covC,
                          "AbC"=AbC,
                          "RMSE"=rmse,
                          "Q"=metares$Q, "tau2"=metares$tau2, "I2"=metares$I2, "H"=metares$H,
                          "anyError"= WarningError[k,2]*1, "ErrorMessage"= WarningError[k,3],
                          "anyWarning"= WarningError[k,4]*1,"WarningMessage"= WarningError[k,5])
      
      ppp[k]<-plogis(metares$TE.random)
      reslist<-list(res=res, pred=ppp,WarningError=WarningError)
    }
  }
  return(reslist)
  
}


###############################################
metareg_logit_tau<-function(meta1, co){  
  #### load true C-index
  load("CIndex_True_Integral_sc1.1_beta0.5_mu0_sd0.5_gumbelLoc0_gumbelScale0.2.rda", verbose=T)
  Ctrue<-C
  Ctrue<-Ctrue[Ctrue$tau>=min(meta1$tau) & Ctrue$tau<=max(meta1$tau),]
  
  #### sort meta-data, calculate inverse variance
  meta1<-meta1[order(meta1$tau),]
  meta1<-helper_inverseVar(meta1)
  temp<-meta1
  meta1<-meta1[!is.na(meta1$SElogitcindex),]  
  
  #### initialize result data frames:
  nr<-rep(NA, nrow(temp))
  ppp<-data.frame("pred_DL_logit_tau"=nr, "se_DL_logit_tau"=nr, 
                  "ci.lb_DL_logit_tau"=nr, "ci.ub_DL_logit_tau"=nr,
                  "pi.lb_DL_logit_tau"=nr, "pi.ub_DL_logit_tau"=nr,
                  "pred_REML_logit_tau"=nr, "se_REML_logit_tau"=nr,
                  "ci.lb_REML_logit_tau"=nr, "ci.ub_REML_logit_tau"=nr,
                  "pi.lb_REML_logit_tau"=nr, "pi.ub_REML_logit_tau"=nr)
  WarningError<-data.frame("var_method"=c(NA,NA),"anyError"=c(NA, NA), "ErrorMessage"=c(NA, NA), 
                           "anyWarning"=c(NA, NA), "WarningMessage"=c(NA, NA))
  res<-rbind(returnNA(name="raw", desc="error", meta1=meta1),
             returnNA(name="raw", desc="error", meta1=meta1))
  
  p1<-1
  methods<-c("DL", "REML")
  for(k in seq_along(methods)){
    #### run meta-regression model
    data <- try(metagen(TE = logitcindex,
                        seTE = SElogitcindex,  studlab = id,
                        data = meta1, sm = "RMD", #SMD
                        comb.fixed = TRUE, comb.random = TRUE,
                        title = "C-Index", method.tau=methods[k]))
    chstr<-paste0("metares<-metareg(x = data, formula =~  tau, method='",methods[k],"', hakn = TRUE)")
    tempres<-withCallingHandlers_LW(chstr = chstr, data=data)
    #### store warnings and errors
    WarningError[k,]<-c(methods[k],tempres$WarningError)
    #### meta-regression model
    metaresRE<-tempres$metares
    
    if(WarningError$anyError[k]){
      res[k,]<-returnNA(name="logit_tau", desc="failed", meta1=meta1)
      res[k,c("anyError",  "ErrorMessage", "anyWarning","WarningMessage")]<-  WarningError[k,2:5]
      reslist<-list(res=res, pred=ppp, WarningError=WarningError)
    }else{
      #### prediction and confidence intervals at 80perc of max. FU time
      newdat<-c("tau"=co)
      pred_80<-as.data.frame(predict(object = metaresRE,newmods = newdat))
      cilower<-plogis(pred_80$ci.lb)  
      ciupper<-plogis(pred_80$ci.ub)  
      
      #### calculate area between prediction and true C-index
      Ctrue$pred<-plogis(predict(object = metaresRE,newmods = as.matrix(Ctrue$tau))$pred)
      AbC<-trapz(x = Ctrue$tau,y = abs(Ctrue$pred-Ctrue$cindex)) / (max(Ctrue$tau)-min(Ctrue$tau))
      
      ####RMSE
      rmse<-sqrt(mean((meta1$cindex-plogis(predict(metaresRE)$pred))^2))
      
      #### Coverage
      trueCindex<-C$cindex[C$tau==co]
      covC<-1*((trueCindex>cilower) & (trueCindex<ciupper))
      
      res[k,]<- data.frame("name"="logit_tau",
                           "description" = "",
                           "var_method"=methods[k],
                           "nstudies"=nrow(meta1),
                           "pred_80"=plogis(pred_80$pred),
                           "se_80"=pred_80$se,
                           "CI_80_lower"=cilower,
                           "CI_80_upper"=ciupper,
                           "inCI"=covC,
                           "AbC"=AbC,
                           "RMSE"=rmse,
                           "Q"=metaresRE$QE, 
                           "tau2"=metaresRE$tau2, 
                           "I2"=metaresRE$I2/100, 
                           "H"=metaresRE$H,
                           "anyError"= WarningError[k,2]*1, "ErrorMessage"= WarningError[k,3],
                           "anyWarning"= WarningError[k,4]*1,"WarningMessage"= WarningError[k,5])
      
      ppp1<-as.data.frame(matrix(NA, nrow=nrow(temp), ncol=6))
      ppp1[!is.na(temp$SElogitcindex) & !is.na(temp$logitcindex),]<-as.data.frame(predict(metaresRE))[,1:6]
      names(ppp1)<-paste0(names(as.data.frame(predict(metaresRE))[,1:6]), "_", methods[k])
      ppp[,p1:(p1+5)]<-ppp1  
      ppp[,p1]<-plogis(ppp[,p1])
      reslist<-list(res=res, pred=ppp, WarningError=WarningError)
      
    }
    p1<-p1+6
  }
  
  return(reslist)
  
}


################################################
metamean_asin<-function(meta1, meta1P, co){  
  #### load true C-index
  load("CIndex_True_Integral_sc1.1_beta0.5_mu0_sd0.5_gumbelLoc0_gumbelScale0.2.rda", verbose=T)
  Ctrue<-C
  Ctrue<-Ctrue[Ctrue$tau>=min(meta1P$tau) & Ctrue$tau<=max(meta1P$tau),]
  
  #### sort meta-data
  meta1<-meta1[order(meta1$tau),]
  meta1<-meta1[!is.na(meta1$SEasincindex),]
  
  #### initialize result data frames:
  data<-data.frame("n"=meta1$n, "cindex"=meta1$asincindex, "sdcindex" = meta1$SEasincindex*sqrt(meta1$n))
  ppp<-c("DL_asin"=NA, "REML_asin"=NA)
  WarningError<-data.frame("var_method"=c(NA,NA),"anyError"=c(NA, NA), "ErrorMessage"=c(NA, NA), 
                           "anyWarning"=c(NA, NA), "WarningMessage"=c(NA, NA))
  res<-rbind(returnNA(name="raw", desc="error", meta1=meta1),
             returnNA(name="raw", desc="error", meta1=meta1))
  
  methods<-c("DL", "REML")
  for(k in seq_along(methods)){
    #### run meta-regression model
    chstr<-paste0("metares<-metamean(n = n,  mean=cindex, sd=sdcindex,data=data, method.tau='",methods[k],"', hakn = TRUE)")
    tempres<-withCallingHandlers_LW(chstr = chstr,data=data)
    #### meta-regression model
    metares<-tempres$metares
    #### store warnings and errors
    WarningError[k,]<-c(methods[k],tempres$WarningError)
    
    if(WarningError$anyError[k]){
      res[k,]<-returnNA(name="asin", desc="failed", meta1=meta1)
      res[k,c("anyError",  "ErrorMessage", "anyWarning","WarningMessage")]<-  WarningError[k,2:5]
      reslist<-list(res=res, pred=ppp, WarningError=WarningError)
    }else{
      #### estimate and confidence intervals
      Ctrue$pred<-sin(metares$TE.random)^2
      cilower<-sin(metares$lower.random)^2
      ciupper<-sin(metares$upper.random)^2
      #### calculate area between prediction and true C-index
      AbC<-trapz(x = Ctrue$tau,y = abs(Ctrue$pred-Ctrue$cindex)) / (max(Ctrue$tau)-min(Ctrue$tau))
      
      #### RMSE
      rmse<-sqrt(mean((meta1P$cindex-sin(metares$TE.random)^2)^2))
      #### coverage
      trueCindex<-C$cindex[C$tau==co]
      covC<-1*((trueCindex>cilower) & (trueCindex<ciupper))
      
      res[k,]<-data.frame("name"="asin",
                          "description" ="",
                          "var_method"=methods[k],
                          "nstudies"=nrow(meta1),
                          "pred_80"=sin(metares$TE.random)^2,
                          "se_80"=metares$seTE.random,
                          "CI_80_lower"=cilower,
                          "CI_80_upper"=ciupper,
                          "inCI"=covC,
                          "AbC"=AbC,
                          "RMSE"=rmse,
                          "Q"=metares$Q, "tau2"=metares$tau2, "I2"=metares$I2, "H"=metares$H,
                          "anyError"= WarningError[k,2]*1, "ErrorMessage"= WarningError[k,3],
                          "anyWarning"= WarningError[k,4]*1,"WarningMessage"= WarningError[k,5])
      
      ppp[k]<-sin(metares$TE.random)^2
      
      reslist<-list(res=res, pred=ppp, WarningError=WarningError)
    }
  }
  
  return(reslist)
  
}


###############################################

metareg_asin_tau<-function(meta1, co){  
  #### load true C-index
  load("CIndex_True_Integral_sc1.1_beta0.5_mu0_sd0.5_gumbelLoc0_gumbelScale0.2.rda", verbose=T)
  Ctrue<-C
  Ctrue<-Ctrue[Ctrue$tau>=min(meta1$tau) & Ctrue$tau<=max(meta1$tau),]

  #### sort meta-data
  meta1<-meta1[order(meta1$tau),]
  meta1<-helper_inverseVar(meta1)
  temp<-meta1
  meta1<-meta1[!is.na(meta1$SEasincindex),]  
  
  #### initialize result data frames:
  nr<-rep(NA, nrow(temp))
  ppp<-data.frame("pred_DL_asin_tau"=nr, "se_DL_asin_tau"=nr, 
                  "ci.lb_DL_asin_tau"=nr, "ci.ub_DL_asin_tau"=nr,
                  "pi.lb_DL_asin_tau"=nr, "pi.ub_DL_asin_tau"=nr,
                  "pred_REML_asin_tau"=nr, "se_REML_asin_tau"=nr,
                  "ci.lb_REML_asin_tau"=nr, "ci.ub_REML_asin_tau"=nr,
                  "pi.lb_REML_asin_tau"=nr, "pi.ub_REML_asin_tau"=nr)
  WarningError<-data.frame("var_method"=c(NA,NA),"anyError"=c(NA, NA), "ErrorMessage"=c(NA, NA), 
                           "anyWarning"=c(NA, NA), "WarningMessage"=c(NA, NA))
  res<-rbind(returnNA(name="raw", desc="error", meta1=meta1),
             returnNA(name="raw", desc="error", meta1=meta1))
  
  methods<-c("DL", "REML")
  p1<-1
  for(k in seq_along(methods)){
    #### run meta-regression model
    data <- try(metagen(TE = asincindex,
                        seTE = SEasincindex,  studlab = id,
                        data = meta1, sm = "RMD", #SMD
                        comb.fixed = TRUE, comb.random = TRUE,
                        title = "C-Index", method.tau=methods[k]))
    chstr<-paste0("metares<-metareg(x = data, formula =~  tau, method='",methods[k],"', hakn = TRUE)")
    tempres<-withCallingHandlers_LW(chstr = chstr,data=data)
    #### store warnings and errors
    WarningError[k,]<-c(methods[k],tempres$WarningError)
    #### meta-regression model
    metaresRE<-tempres$metares
    
    
    if(WarningError$anyError[k]){
      res[k,]<-returnNA(name="asin_tau", desc="failed", meta1=meta1)
      res[k,c("anyError",  "ErrorMessage", "anyWarning","WarningMessage")]<-  WarningError[k,2:5]
      reslist<-list(res=res, pred=ppp, WarningError=WarningError)
    }else{
      
      #### estimate and confidence intervals and 80perc. of max. FU time
      newdat<-c("tau"=co)
      pred_80<-as.data.frame(predict(object = metaresRE,newmods = newdat))
     
      cilower<- sin(pred_80$ci.lb)^2 #pred_80[1]-qt(0.975, df=nrow(meta1)-2)*pred_80[2]
      ciupper<- sin(pred_80$ci.ub)^2 #pred_80[1]+qt(0.975, df=nrow(meta1)-2)*pred_80[2]
      
      #### calculate area between prediction and true C-index
      Ctrue$pred<-sin(predict(object = metaresRE,newmods = as.matrix(Ctrue$tau))$pred)^2
      AbC<-trapz(x = Ctrue$tau,y = abs(Ctrue$pred-Ctrue$cindex)) / (max(Ctrue$tau)-min(Ctrue$tau))
      
      #### RMSE
      rmse<-sqrt(mean((meta1$cindex-sin(predict(metaresRE, level = 0)$pred)^2)^2))
      
      #### Coverage
      trueCindex<-C$cindex[C$tau==co]
      covC<-1*((trueCindex>cilower) & (trueCindex<ciupper))
      
      res[k,]<- data.frame("name"="asin_tau",
                           "description" = "",
                           "var_method"=methods[k],
                           "nstudies"=nrow(meta1),
                           "pred_80"=sin(pred_80$pred)^2,
                           "se_80"=pred_80$se,
                           "CI_80_lower"=cilower,
                           "CI_80_upper"=ciupper,
                           "inCI"=covC,
                           "AbC"=AbC,
                           "RMSE"=rmse,
                           "Q"=metaresRE$QE, 
                           "tau2"=metaresRE$tau2, 
                           "I2"=metaresRE$I2/100, 
                           "H"=metaresRE$H,
                           "anyError"= WarningError[k,2]*1, "ErrorMessage"= WarningError[k,3],
                           "anyWarning"= WarningError[k,4]*1,"WarningMessage"= WarningError[k,5])
      
      ppp1<-as.data.frame(matrix(NA, nrow=nrow(temp), ncol=6))
      ppp1[!is.na(temp$SEasincindex) & !is.na(temp$asincindex),]<-as.data.frame(predict(metaresRE))[,1:6]
      names(ppp1)<-paste0(names(as.data.frame(predict(metaresRE))[,1:6]), "_", methods[k])
      ppp[,p1:(p1+5)]<-ppp1  
      ppp[,p1]<-sin(ppp[,p1])^2
      
      reslist<-list(res=res, pred=ppp, WarningError=WarningError)
      
    }
    p1<-p1+6
  }
  
  return(reslist)
  
  
}


###############################################
metareg_tau_spline<-function(meta1, cindex="cindex", SEcindex="SEcindex", co){  
  #### load true C-index
  load("CIndex_True_Integral_sc1.1_beta0.5_mu0_sd0.5_gumbelLoc0_gumbelScale0.2.rda", verbose=T)
  #### temporary data frames meta2 and meta3
  meta2<-meta3<-C
  wwmeta2<-which(meta2$tau==co)
  meta2<-meta2[meta2$tau>=min(meta1$tau) & meta2$tau<=max(meta1$tau),]
  
  #### sort meta data
  meta1<-meta1[order(meta1$tau),]
  
  #### initialize result data frames:
  nr<-rep(NA, nrow(meta1))
  ppp<-data.frame("pred_DL_rcs_tau"=nr, "se_DL_rcs_tau"=nr, 
                  "ci.lb_DL_rcs_tau"=nr, "ci.ub_DL_rcs_tau"=nr,
                  "pi.lb_DL_rcs_tau"=nr, "pi.ub_DL_rcs_tau"=nr,
                  "pred_REML_rcs_tau"=nr, "se_REML_rcs_tau"=nr,
                  "ci.lb_REML_rcs_tau"=nr, "ci.ub_REML_rcs_tau"=nr,
                  "pi.lb_REML_rcs_tau"=nr, "pi.ub_REML_rcs_tau"=nr)
  WarningError<-data.frame("var_method"=c(NA,NA),"anyError"=c(NA, NA), "ErrorMessage"=c(NA, NA), 
                           "anyWarning"=c(NA, NA), "WarningMessage"=c(NA, NA))
  res<-rbind(returnNA(name="raw", desc="error", meta1=meta1),
             returnNA(name="raw", desc="error", meta1=meta1))
  
  methods<-c("DL", "REML")
  p1<-1
  for(k in seq_along(methods)){
    
    nknots<-4
    if(nrow(meta1)<30) nknots<-3
    
    #### run meta-regression model
     knotsrma<-rcs(meta1$tau, nknots) 
    chstr<-paste0("metares<-rma(yi = ",cindex,", vi = ",
                  SEcindex,"^2, mods=~rcs(tau, parms=c(", 
                  paste0(attr(knotsrma, "parms"), collapse = ", "),
                  ")), data=data, method='",methods[k],"', test='knha')")
    tempres<-withCallingHandlers_LW(chstr = chstr,data=meta1)
    #### store warnings and errors
    WarningError[k,]<-c(methods[k],tempres$WarningError)
    #### meta-regression model
    metaresRE<-tempres$metares

    if(WarningError$anyError[k]){
      res[k,]<-returnNA(name=paste0("rmaspline_", cindex), desc="failed", meta1=meta1)
      res$var_method[k]<-methods[k]
      res[k,c("anyError",  "ErrorMessage", "anyWarning","WarningMessage")]<-  WarningError[k,2:5]
      reslist<-list(res=res, pred=ppp, WarningError=WarningError)
    }else{
      
      #### heterogeneity measures
      tau2est<-metaresRE$tau2
      I2est<-metaresRE$I2
      Hest<-metaresRE$H2
      CQ<-metaresRE$QE
      
      ### predictions at 80perc of max. FU time
      newdat<-data.frame("tau"=c(meta1$tau, co))
      X <- model.matrix(~ rcs(tau, attr(knotsrma, "parms")), data=newdat)[,-1]
      colnames(X)<-names(coef(metaresRE))[-1]
      pred_80<-as.data.frame(predict(object = metaresRE,newmods = X))
      pred_80<-as.numeric(pred_80[nrow(pred_80),1:2])

      X <- model.matrix(~ rcs(tau, attr(knotsrma, "parms")), data=meta2)[,-1]
      colnames(X)<-names(coef(metaresRE))[-1]
      meta2$pred<-predict(object = metaresRE,newmods = X)$pred
      meta2$se<-predict(object = metaresRE,newmods = X)$se
      meta2$ci.lb<-predict(object = metaresRE,newmods = X)$ci.lb
      meta2$ci.ub<-predict(object = metaresRE,newmods = X)$ci.ub
      
      X <- model.matrix(~ rcs(tau, attr(knotsrma, "parms")), data=meta3)[,-1]
      colnames(X)<-names(coef(metaresRE))[-1]
      meta3$ci.lb<-predict(object = metaresRE,newmods = X)$ci.lb
      meta3$ci.ub<-predict(object = metaresRE,newmods = X)$ci.ub
      
      #### predictions on training data
      predmeta1<-predict(metaresRE)$pred    
      
      #### confidence intervals 
      trueCindex<-C$cindex[C$tau==co]
      cilower<-meta3$ci.lb[meta3$tau==co]
      ciupper<-meta3$ci.ub[meta3$tau==co]

      
      if(cindex=="logitcindex"){
        pred_80[1]<-plogis(pred_80[1]); 
        meta2$pred<-plogis(meta2$pred)
        predmeta1<-plogis(predmeta1)
        cilower<-plogis(cilower)
        ciupper<-plogis(ciupper)
        }

      if(cindex=="asincindex"){
        pred_80[1]<-sin(pred_80[1])^2;
        meta2$pred<-sin(meta2$pred)^2
        predmeta1<-sin(predmeta1)^2
        cilower<-sin(cilower)^2
        ciupper<-sin(ciupper)^2
      }
      
      #### coverage
      covC<-1*((trueCindex>cilower)& (trueCindex<ciupper))
      
      #### calculate area between prediction and true C-index
      AbC<-trapz(x = meta2$tau,y = abs(meta2$pred-meta2$cindex)) / (max(meta2$tau)-min(meta2$tau))
      rmse<-sqrt(mean((meta1$cindex-predmeta1)^2))
      
      res[k,]<-data.frame("name"=paste0("rmaspline_", cindex),
                          "description" = "",
                          "var_method"=methods[k],
                          "nstudies"=nrow(meta1),
                          "pred_80"=pred_80[1],
                          "se_80"=pred_80[2],
                          "CI_80_lower"=cilower,
                          "CI_80_upper"=ciupper,
                          "inCI"=covC,
                          "AbC"=AbC,
                          "RMSE"=rmse,
                          "Q"=CQ,
                          "tau2"=tau2est,
                          "I2"=I2est,
                          "H"=Hest,
                          "anyError"= WarningError[k,2]*1, "ErrorMessage"= WarningError[k,3],
                          "anyWarning"= WarningError[k,4]*1,"WarningMessage"= WarningError[k,5])
      
      ppp1<-as.data.frame(matrix(NA, nrow=nrow(meta1), ncol=6))
      ppp1[!is.na(meta1[,cindex]) & !is.na(meta1[,SEcindex]),]<-as.data.frame(predict(metaresRE, level=0))[,1:6]
      names(ppp1)<-paste0(names(as.data.frame(predict(metaresRE))[,1:6]), "_", methods[k])
      ppp[,p1:(p1+5)]<-ppp1  
      
      if(cindex=="logitcindex"){ppp[,p1]<-plogis(ppp[,p1])}
      if(cindex=="asincindex"){ppp[,p1]<-sin(ppp[,p1])^2}
      
      reslist<-list(res=res, pred=ppp, WarningError=WarningError)
    }
    p1<-p1+6
  }
  
  return(reslist)
  
}




###############################################
metamean_mfp<-function(meta1, cindex="cindex", SEcindex="SEcindex", co){  
  #### load true C-index
  load("CIndex_True_Integral_sc1.1_beta0.5_mu0_sd0.5_gumbelLoc0_gumbelScale0.2.rda", verbose=T)
  
  
  #### sort meta-data
  meta1<-meta1[order(meta1$tau),]
  #### calculate fractional polynomials
  meta1$var1<-meta1$tau^{-0.5}
  meta1$var2<-meta1$tau^{0.5}
  
  #### temporary data frames meta2 and meta3
  meta2<-meta3<-C
  meta2<-meta2[meta2$tau>=min(meta1$tau) & meta2$tau<=max(meta1$tau),]
  wwmeta2<-which(meta2$tau==co)
  meta2$var1<-meta2$tau^{-0.5}
  meta2$var2<-meta2$tau^{0.5}
  meta3$var1<-meta3$tau^{-0.5}
  meta3$var2<-meta3$tau^{0.5}
  newdata<-as.matrix(meta2[,c(3,4)])
  newdata3<-as.matrix(meta3[,c(3,4)])

  
  #### initialize result data frames:
  nr<-rep(NA, nrow(meta1))
  ppp<-data.frame("pred_DL_mfp_tau"=nr, "se_DL_mfp_tau"=nr, 
                  "ci.lb_DL_mfp_tau"=nr, "ci.ub_DL_mfp_tau"=nr,
                  "pi.lb_DL_mfp_tau"=nr, "pi.ub_DL_mfp_tau"=nr,
                  "pred_REML_mfp_tau"=nr, "se_REML_mfp_tau"=nr,
                  "ci.lb_REML_mfp_tau"=nr, "ci.ub_REML_mfp_tau"=nr,
                  "pi.lb_REML_mfp_tau"=nr, "pi.ub_REML_mfp_tau"=nr)
  WarningError<-data.frame("var_method"=c(NA,NA),"anyError"=c(NA, NA), "ErrorMessage"=c(NA, NA), 
                           "anyWarning"=c(NA, NA), "WarningMessage"=c(NA, NA))
  res<-rbind(returnNA(name="raw", desc="error", meta1=meta1),
             returnNA(name="raw", desc="error", meta1=meta1))
  

  
  methods<-c("DL", "REML")
  p1<-1
  for(k in seq_along(methods)){
    #### run meta-regression model
    chstr<-paste0("metares<-rma(yi = ",cindex,", vi = ",
                  SEcindex,"^2, mods=~var1+var2, data=data, method='",
                  methods[k],"', test='knha')")
    tempres<-withCallingHandlers_LW(chstr = chstr,data=meta1)
    #### store warnings and errors
    WarningError[k,]<-c(methods[k],tempres$WarningError)
    #### meta-regression model
    metaresRE<-tempres$metares
    
    if(WarningError$anyError[k]){
      res[k,]<-returnNA(name=paste0("fracpoly_",cindex), desc="failed", meta1=meta1)
      res$var_method[k]<-methods[k]
      res[k,c("anyError",  "ErrorMessage", "anyWarning","WarningMessage")]<-  WarningError[k,2:5]
      reslist<-list(res=res, pred=ppp, WarningError=WarningError)
    }else{
      
      #### heterogeneity measures
      tau2est<-metaresRE$tau2
      I2est<-metaresRE$I2
      Hest<-metaresRE$H2
      CQ<-metaresRE$QE
      
      #### predictions at 80perc. of max. follow-up times
      pred_newmod<-as.data.frame(predict(object = metaresRE,newmods = newdata))
      pred_80<-as.numeric(pred_newmod[wwmeta2,1:2])
      
      meta2$pred<-pred_newmod$pred
      meta2$se<-pred_newmod$se
      meta2$ci.lb<-pred_newmod$ci.lb
      meta2$ci.ub<-pred_newmod$ci.ub
      
      pred_newmod3<-as.data.frame(predict(object = metaresRE,newmods = newdata3))
      meta3$ci.lb<-pred_newmod3$ci.lb
      meta3$ci.ub<-pred_newmod3$ci.ub
      
      #### confidence intervals
      cilower<-meta3$ci.lb[meta3$tau==co]
      ciupper<-meta3$ci.ub[meta3$tau==co]

      predmeta1<-predict(metaresRE)$pred    
      trueCindex<-C$cindex[C$tau==co]
      
      
      if(cindex=="logitcindex"){
        cilower<-plogis(cilower)
        ciupper<-plogis(ciupper)
        pred_80[1]<-plogis(pred_80[1]); 
        meta2$pred<-plogis(meta2$pred)
        predmeta1<-plogis(predmeta1)

      }

      if(cindex=="asincindex"){
        cilower<-sin(cilower)^2
        ciupper<-sin(ciupper)^2
        pred_80[1]<-sin(pred_80[1])^2;
        meta2$pred<-sin(meta2$pred)^2
        predmeta1<-sin(predmeta1)^2
      }

      #### coverage
      covC<-1*((trueCindex>cilower)& (trueCindex<ciupper))
      #### calculate area between prediction and true C-index
      AbC<-trapz(x = meta2$tau,y = abs(meta2$pred-meta2$cindex)) / (max(meta2$tau)-min(meta2$tau))
      ####RMSE
      rmse<-sqrt(mean((meta1$cindex-predmeta1)^2))
      
      res[k,]<-data.frame("name"=paste0("fracpoly_", cindex),
                          "description" = "",
                          "var_method"=methods[k],
                          "nstudies"=nrow(meta1),
                          "pred_80"=pred_80[1],
                          "se_80"=pred_80[2],
                          "CI_80_lower"=cilower,
                          "CI_80_upper"=ciupper,
                          "inCI"=covC,
                          "AbC"=AbC,
                          "RMSE"=rmse,
                          "Q"=CQ,
                          "tau2"=tau2est,
                          "I2"=I2est,
                          "H"=Hest,
                          "anyError"= WarningError[k,2]*1, "ErrorMessage"= WarningError[k,3],
                          "anyWarning"= WarningError[k,4]*1,"WarningMessage"= WarningError[k,5])
      
      ppp1<-as.data.frame(matrix(NA, nrow=nrow(meta1), ncol=6))
      ppp1[!is.na(meta1[,cindex]) & !is.na(meta1[,SEcindex]),]<-as.data.frame(predict(metaresRE))[,1:6]
      names(ppp1)<-paste0(names(as.data.frame(predict(metaresRE))[,1:6]), "_", methods[k])
      ppp[,p1:(p1+5)]<-ppp1  
      
      if(cindex=="logitcindex"){ppp[,p1]<-plogis(ppp[,p1])}
      if(cindex=="asincindex"){ppp[,p1]<-sin(ppp[,p1])^2}
      
      reslist<-list(res=res, pred=ppp, WarningError=WarningError)
    }
    p1<-p1+6
  }
  
  return(reslist)
  
}

########################################################################################
metamean_expdecay_lw3_R0free<-function(meta1,  cindex="cindex", invVar="invVar", co){  
  #### load true C-index
  load("CIndex_True_Integral_sc1.1_beta0.5_mu0_sd0.5_gumbelLoc0_gumbelScale0.2.rda", verbose=T)
  Ctrue<-C
  Ctrue<-Ctrue[Ctrue$tau>=min(meta1$tau) & Ctrue$tau<=max(meta1$tau),]
  
  #### sort meta-data
  meta1<-meta1[order(meta1$tau),]
  meta1<-helper_inverseVar(meta1)
  temp<-meta1
  meta1<-meta1[!is.na(meta1[,invVar]),]

  #################################################################################
  #### run meta-regression model

  data <- groupedData(cindex ~ tau | id, meta1)
  formula_nlme<-paste0("metares<-nlme(",cindex," ~ SSasymp(tau, Asymp, R0, lrc),
                    fixed = R0 +Asymp + lrc~ 1,
                    random =  Asymp~ 1|id, 
                    weights = varFixed(~ 1/",invVar,"),
                    data = data, start = c(R0=1, Asymp = 0.75, lrc=1),
                    method='REML',
                    control = nlmeControl(sigma = 1, returnObject = T))")
  
  chstr<-gsub("\n","",formula_nlme)
  tempres<-withCallingHandlers_LW(chstr = chstr, data=data)
  #### store warnings and errors
  WarningError<-tempres$WarningError
  #### meta-regression model
  metaresRE<-tempres$metares
  
  bestfit<-tempres
  temp_WarningError<-WarningError
  temp_metaresRE<-metaresRE
  
  res<-rbind(returnNA(name="raw", desc="error", meta1=meta1))
  
  ### try to find better solution if warning or error by
  #### grid search of starting values:
  if(WarningError$anyError | WarningError$anyWarning){
      
      start_nlme <-expand.grid( R0=seq(0.2, 2, 0.4), Asymp=seq(0.6,0.9, 0.06), lrc=seq(0.5,3,0.4))
      start_nlme<-start_nlme[sample(1:nrow(start_nlme)),]
      
      boolstop<-convbestfit<-FALSE
      k<-1
      maxit<-50 # maximum iteration number
      noconv<-vector()
  
      #### while warning or error, or maximum number of iteration is not reached
      while(!boolstop & k<= nrow(start_nlme)){
        #### starting values
        a<-start_nlme[k,2]
        lrc<-start_nlme[k,3]
        r0<-start_nlme[k,1]
        
        print(paste0("R0=", r0,  "; lrc: ", lrc, "; Asymptote: ", a))
        #### initalize model
        formula_nlme<-paste0("metares<-nlme(",cindex," ~ SSasymp(tau, Asymp, R0, lrc),
                        fixed = R0 +Asymp + lrc~ 1,
                        random =  Asymp~ 1|id,
                        weights = varFixed(~ 1/",invVar,"),
                        data = data, start = c(R0=",r0,", Asymp = ",a ,", lrc=",lrc,"),
                        method='REML',
                        control = nlmeControl(sigma = 1, maxIter = ",maxit,", returnObject = T))")
        
        chstr<-gsub("\n","",formula_nlme)
        
        tryCatch({
          #### run model
        setTimeLimit(60) ## set time limit for function
        tempres<-withCallingHandlers_LW(chstr = chstr, data=data)

        ## if no error, no warning -> stop
        if(!tempres$WarningError$anyError & !tempres$WarningError$anyWarning){ 
          bestfit<-tempres
          temp_WarningError<-tempres$WarningError
          boolstop<-TRUE
          #### if warning but no error -> continue
        } else if (!tempres$WarningError$anyError & tempres$WarningError$anyWarning &
                   temp_WarningError$anyError){
          bestfit<-tempres
          temp_WarningError<-tempres$WarningError
          #### if error -> continue
        } else if (!tempres$WarningError$anyError & tempres$WarningError$anyWarning &  ## no error but warning now
                   grepl("reached without convergence", temp_WarningError$WarningMessage) & ## no convergence in last bestfit
                   !grepl("reached without convergence", tempres$WarningError$WarningMessage)){##no convergence issue in this run
          bestfit<-tempres
          temp_WarningError<-tempres$WarningError
        }
        
        if(grepl("reached without convergence", tempres$WarningError$WarningMessage)){
          noconv<-c(noconv, k)
        }
        }, error = function(e) NULL)
        k<-k+1
      }
      
      if(!grepl("reached without convergence", tempres$WarningError$WarningMessage)){
        convbestfit<-TRUE
      }
  
    }
  
  WarningError<-data.frame("var_method"="REML", bestfit$WarningError)
  metaresRE<-bestfit$metares
  
  #### results:
  if(!(WarningError$anyError)){
    coefnlsRE<-tau2est<-I2est<-Hest<-pFE<-pRE<-NA
    vc <- VarCorr(metaresRE)
    suppressWarnings(storage.mode(vc) <- "numeric")
    tau2est<-vc[1,"StdDev"]^2
    
    #### predictions for 80perc of max. FU
    newdat<-data.frame("tau"=co)
    pred_80<-as.numeric(as.data.frame(predict(object = metaresRE,level=0,newdata = newdat)))
    
    Ctrue$pred<-as.data.frame(predict(object = metaresRE,newdata = Ctrue, level=0))[,1] ## for Area btw. curves
    predmeta1<-predict(metaresRE, level=0)
    
    if(cindex=="logitcindex"){
      pred_80<-plogis(pred_80); 
      Ctrue$pred<-plogis(Ctrue$pred)
      predmeta1<-plogis(predmeta1)
    }
    if(cindex=="logcindex"){
      pred_80<-exp(pred_80); 
      Ctrue$pred<-exp(Ctrue$pred)
      predmeta1<-exp(predmeta1)
    }
    if(cindex=="asincindex"){
      pred_80<-sin(pred_80)^2; 
      Ctrue$pred<-sin(Ctrue$pred)^2
      predmeta1<-sin(predmeta1)^2
    }
    
    #### calculate area between prediction and true C-index
    AbC<-trapz(x = Ctrue$tau,y = abs(Ctrue$pred-Ctrue$cindex)) / (max(Ctrue$tau)-min(Ctrue$tau))
    ####RMSE
    rmse<-sqrt(mean((meta1$cindex-predmeta1)^2))
    
    cilower<-ciupper<-covC<-NA
    
    res<-data.frame("name"=paste0("expdecay_",cindex),
                    "description" = "",
                    "var_method"="REML",
                    "nstudies"=nrow(meta1), 
                    "pred_80"=pred_80,
                    "se_80"=NA,
                    "CI_80_lower"=cilower,
                    "CI_80_upper"=ciupper,
                    "inCI"=covC,
                    "AbC"=AbC,
                    "RMSE"=rmse,
                    "Q"=NA, "tau2"=tau2est, "I2"=NA, "H"=NA,
                    "anyError"= WarningError[2]*1, "ErrorMessage"= WarningError[3],
                    "anyWarning"= WarningError[4]*1,"WarningMessage"= WarningError[5])
    
    
    
    paramsRE<-fixef(metaresRE)
    
    ppp<-as.data.frame(matrix(NA, nrow=nrow(temp), ncol=1))
    ppp[!is.na(temp[,cindex]) & !is.na(temp[,invVar]),]<-as.data.frame(predict(metaresRE, level=0))
    names(ppp)[1]<-"pred_exponentialdecay_REML"
    ppp$pred_exponentialdecay_asymp_REML<-paramsRE[1]
    
    if(cindex=="logitcindex"){
      ppp$pred_exponentialdecay_REML<-plogis(ppp$pred_exponentialdecay_REML);
      ppp$pred_exponentialdecay_asymp_REML<-plogis(ppp$pred_exponentialdecay_asymp_REML);
    }
    if(cindex=="logcindex"){
      ppp$pred_exponentialdecay_REML<-exp(ppp$pred_exponentialdecay_REML);
      ppp$pred_exponentialdecay_asymp_REML<-exp(ppp$pred_exponentialdecay_asymp_REML);
    }
    if(cindex=="asincindex"){
      ppp$pred_exponentialdecay_REML<-sin(ppp$pred_exponentialdecay_REML)^2;
      ppp$pred_exponentialdecay_asymp_REML<-sin(ppp$pred_exponentialdecay_asymp_REML)^2;
    }
    
    reslist<-list(res=res, pred=ppp, WarningError=WarningError)
    
  }else{
    ppp<-as.data.frame(matrix(NA, nrow=nrow(temp), ncol=2))
    names(ppp)<-c("pred_exponentialdecay_REML", "pred_exponentialdecay_asymp_REML")
    
    res[1,]<-returnNA(name=paste0("expdecay_", cindex), desc="failed", meta1=meta1)
    res$var_method[1]<-"REML"
    res[1,c("anyError",  "ErrorMessage", "anyWarning","WarningMessage")]<-  WarningError[2:5]
    reslist<-list(res=res, pred=ppp, WarningError=WarningError)
  }
  
  return(reslist)
  
}
