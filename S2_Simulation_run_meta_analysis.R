library("batchtools")

## function to  run all meta-analysis and -regression models
###################################################################################
####

#### file: file with meta-analysis data (generated from XXX)
#### nstudies_per_meta: number of studies included in meta-analysis
#### ww: simulation run

file<-"simulation/Simulation_generate_meta_data_1.rda"
cind_meta<-function(nstudies_per_meta, ww, file){
  source("S0_functions.R")
  
  #### load data
  load(file)
  studies<-output$studies
  setting<-as.vector(output$setting[ww,])
  names(setting)<-colnames(output$setting)
  
  #### 80percent of follow-up time
  co80<-round(0.8*as.numeric(setting["taumax"]), 2) ## Auswertezeitpunkt
  
  #### sort data
  meta1<-studies[[ww]]
  meta1<-meta1[1:nstudies_per_meta,]
  meta1<-meta1[!is.na(meta1$cindex) & meta1$SEcindex>0 & !is.na(meta1$SElogitcindex),]
  meta1<-meta1[order(meta1$tau),]
  
  #### data frames containing only 50 or 30 percent of the data:
  meta_last50perc<-meta1[round(nrow(meta1)/2+1):nrow(meta1),]
  meta_last30perc<-meta1[round(2*nrow(meta1)/3+1):nrow(meta1),]
  
  #### initalize results data frame
  methodnames<-c('raw', 'raw_last50perc', 'raw_last30perc', 
                'logit', 'logit_last50perc', 'logit_last30perc',
                 'asin',  'asin_last50perc', 'asin_last30perc',
                'raw_tau', 'logit_tau', 'asin_tau', 
                'rmaspline_cindex', 'rmaspline_logitcindex', 
                'rmaspline_asincindex', 
                'fracpoly_cindex', 'fracpoly_logitcindex', 
                'fracpoly_asincindex',
                'expdecay_cindex', 'expdecay_logitcindex',
                 'expdecay_asincindex')
  
  results<-data.frame(name=factor(39, levels = methodnames),
                      description=factor(39, levels = c('', 'failed')),
                      var_method=factor(39, levels=c("DL", "REML")), 
                      nstudies=numeric(39), 
                      pred_80=numeric(39), se_80=numeric(39),CI_80_lower=numeric(39),
                      CI_80_upper=numeric(39),
                      inCI=numeric(39),
                      AbC=numeric(39),
                      RMSE=numeric(39),
                      Q=numeric(39), tau2=numeric(39), I2=numeric(39), H=numeric(39),
                      anyError=numeric(39),ErrorMessage=rep(NA,39),
                      anyWarning=numeric(39),WarningMessage=rep(NA,39))
  
  ##########################################################################
  #### time-independent methods:
  ##########################################################################
  print("raw")
  we<-1
  nn<-1
  ####------------------------------------------------------------
  #### untransformed
  model_metamean_raw<-metamean_raw(meta1, meta1P=meta1, co = co80)
  results[nn:(nn+1),]<-model_metamean_raw$res; 
  results$name[nn:(nn+1)]<-"raw"; nn<-nn+2
  meta1$pred_raw_DL<-model_metamean_raw$pred[1]
  meta1$pred_raw_REML<-model_metamean_raw$pred[2]
  ####------------------------------------------------------------
  print("raw")
  #### untransformed - last 50 %
  model_metamean_raw_last50<-metamean_raw(meta1=meta_last50perc, meta1P=meta1, co = co80)
  results[nn:(nn+1),]<-model_metamean_raw_last50$res
  results$name[nn:(nn+1)]<-"raw_last50perc"; nn<-nn+2
  meta1$pred_raw_last50_DL<-model_metamean_raw_last50$pred[1]
  meta1$pred_raw_last50_REML<-model_metamean_raw_last50$pred[2]
  ####------------------------------------------------------------
  #### untransformed - last 30 %
  print("raw")
  model_metamean_raw_last30<-metamean_raw(meta_last30perc, meta1P=meta1, co = co80)
  results[nn:(nn+1),]<-model_metamean_raw_last30$res
  results$name[nn:(nn+1)]<-"raw_last30perc"; nn<-nn+2
  meta1$pred_raw_last30_DL<-model_metamean_raw_last30$pred[1]
  meta1$pred_raw_last30_REML<-model_metamean_raw_last30$pred[2]
  ####################################################
  #### logit-transformed
  print("logit")
  model_metamean_logit<-metamean_logit(meta1, meta1P=meta1, co = co80)
  results[nn:(nn+1),]<-model_metamean_logit$res; 
  results$name[nn:(nn+1)]<-"logit"; nn<-nn+2
  meta1$pred_logit_DL<-model_metamean_logit$pred[1]
  meta1$pred_logit_REML<-model_metamean_logit$pred[2]
  ####------------------------------------------------------------
  #### logit-transformed - last 50%
  print("logit")
  model_metamean_logit_last50<-metamean_logit(meta_last50perc, meta1P=meta1, co = co80)
  results[nn:(nn+1),]<-model_metamean_logit_last50$res
  results$name[nn:(nn+1)]<-"logit_last50perc"; nn<-nn+2
  meta1$pred_logit_last50_DL<-model_metamean_logit_last50$pred[1]
  meta1$pred_logit_last50_REML<-model_metamean_logit_last50$pred[2]
  ####------------------------------------------------------------
  #### logit-transformed - last 30%
  print("logit")
  model_metamean_logit_last30<-metamean_logit(meta_last30perc, meta1P=meta1, co = co80)
  results[nn:(nn+1),]<-model_metamean_logit_last30$res
  results$name[nn:(nn+1)]<-"logit_last30perc"; nn<-nn+2
  meta1$pred_logit_last30_DL<-model_metamean_logit_last30$pred[1]
  meta1$pred_logit_last30_REML<-model_metamean_logit_last30$pred[2]
   ####################################################
  #### asin-transformed 
  print("asin")
  model_metamean_asin<-metamean_asin(meta1, meta1P=meta1, co = co80)
  results[nn:(nn+1),]<-model_metamean_asin$res
  results$name[nn:(nn+1)]<-"asin"; nn<-nn+2
  meta1$pred_asin_DL<-model_metamean_asin$pred[1]
  meta1$pred_asin_REML<-model_metamean_asin$pred[2]
  ####---------------------------------------------
  #### asin-transformed - last 50%
  print("asin")
  model_metamean_asin_last50<-metamean_asin(meta_last50perc, meta1P=meta1, co = co80)
  results[nn:(nn+1),]<-model_metamean_asin_last50$res
  results$name[nn:(nn+1)]<-"asin_last50perc"; nn<-nn+2
  meta1$pred_asin_last50_DL<-model_metamean_asin_last50$pred[1]
  meta1$pred_asin_last50_REML<-model_metamean_asin_last50$pred[2]
  ####---------------------------------------------
  #### asin-transformed  - last 30%
  print("asin")
  model_metamean_asin_last30<-metamean_asin(meta_last30perc, meta1P=meta1, co = co80)
  results[nn:(nn+1),]<-model_metamean_asin_last30$res
  results$name[nn:(nn+1)]<-"asin_last30perc"; nn<-nn+2
  meta1$pred_asin_last30_DL<-model_metamean_asin_last30$pred[1]
  meta1$pred_asin_last30_REML<-model_metamean_asin_last30$pred[2]
  
  ##########################################################################
  #### time dependent meta-regression models
  ##########################################################################
  
  #### linear meta-regression (untransformed C-index)
  print("raw_tau")
  model_metamean_raw_tau<-metareg_raw_tau(meta1, co=co80)
  results[nn:(nn+1),]<-model_metamean_raw_tau$res
  results$name[nn:(nn+1)]<-"raw_tau"; nn<-nn+2
  meta1$pred_raw_tau_DL<-model_metamean_raw_tau$pred$pred_DL_raw_tau
  meta1$se_raw_tau_DL<-model_metamean_raw_tau$pred$se_DL_raw_tau
  meta1$pred_raw_tau_REML<-model_metamean_raw_tau$pred$pred_REML_raw_tau
  meta1$se_raw_tau_REML<-model_metamean_raw_tau$pred$se_REML_raw_tau
  ####---------------------------------------------
  #### linear meta-regression (logit-transformed C-index)
  print("logit_tau")
  model_metareg_logit_tau<-metareg_logit_tau(meta1, co=co80)
  results[nn:(nn+1),]<-model_metareg_logit_tau$res; 
  results$name[nn:(nn+1)]<-"logit_tau"; nn<-nn+2
  meta1$pred_logit_tau_DL[!is.na(meta1$SElogitcindex)]<-model_metareg_logit_tau$pred$pred_DL_logit_tau
  meta1$se_logit_tau_DL[!is.na(meta1$SElogitcindex)]<-model_metareg_logit_tau$pred$se_DL_logit_tau
  meta1$pred_logit_tau_REML[!is.na(meta1$SElogitcindex)]<-model_metareg_logit_tau$pred$pred_REML_logit_tau
  meta1$se_logit_tau_REML[!is.na(meta1$SElogitcindex)]<-model_metareg_logit_tau$pred$se_REML_logit_tau

  ####---------------------------------------------
  #### linear meta-regression (asin-transformed C-index)
  print("asin_tau")
  model_metareg_asin_tau<-metareg_asin_tau(meta1, co=as.numeric(co80))
  results[nn:(nn+1),]<-model_metareg_asin_tau$res; 
  results$name[nn:(nn+1)]<-"asin_tau"; nn<-nn+2
  meta1$pred_asin_tau_DL<-model_metareg_asin_tau$pred$pred_DL_asin_tau
  meta1$se_asin_tau_DL<-model_metareg_asin_tau$pred$se_DL_asin_tau
  meta1$pred_asin_tau_REML<-model_metareg_asin_tau$pred$pred_REML_asin_tau
  meta1$se_asin_tau_REML<-model_metareg_asin_tau$pred$se_REML_asin_tau
  
  ####################################################
  #### spline meta-regression (untransformed C-index)
  print("spline_tau")
  model_metareg_tau_spline<-metareg_tau_spline(meta1, co=co80)
  results[nn:(nn+1),]<-model_metareg_tau_spline$res; 
  results$name[nn:(nn+1)]<-"rmaspline_cindex"; nn<-nn+2
  meta1$pred_rcs_tau_DL<-model_metareg_tau_spline$pred$pred_DL_rcs_tau
  meta1$se_rcs_tau_DL<-model_metareg_tau_spline$pred$se_DL_rcs_tau
  meta1$pred_rcs_tau_REML<-model_metareg_tau_spline$pred$pred_REML_rcs_tau
  meta1$se_rcs_tau_REML<-model_metareg_tau_spline$pred$se_REML_rcs_tau
  ####---------------------------------------------
  #### spline meta-regression  (logit-transformed C-index)
  print("spline_tau, logit")
  model_metareg_tau_spline_logit<-metareg_tau_spline(meta1, cindex = "logitcindex",SEcindex = "SElogitcindex", co=co80)
  results[nn:(nn+1),]<-model_metareg_tau_spline_logit$res; 
  results$name[nn:(nn+1)]<-"rmaspline_logitcindex"; nn<-nn+2
  meta1$pred_rcs_logit_tau_DL<-model_metareg_tau_spline_logit$pred$pred_DL_rcs_tau
  meta1$se_rcs_logit_tau_DL<-model_metareg_tau_spline_logit$pred$se_DL_rcs_tau
  meta1$pred_rcs_logit_tau_REML<-model_metareg_tau_spline_logit$pred$pred_REML_rcs_tau
  meta1$se_rcs_tau_logit_REML<-model_metareg_tau_spline_logit$pred$se_REML_rcs_tau
   ####---------------------------------------------
  #### spline meta-regression (asin-transformed C-index)
  print("spline_tau, asin")
  model_metareg_tau_spline_asin<-metareg_tau_spline(meta1, cindex = "asincindex",SEcindex = "SEasincindex", co=co80)
  results[nn:(nn+1),]<-model_metareg_tau_spline_asin$res; 
  results$name[nn:(nn+1)]<-"rmaspline_asincindex"; nn<-nn+2
  meta1$pred_rcs_asin_tau_DL<-model_metareg_tau_spline_asin$pred$pred_DL_rcs_tau
  meta1$se_rcs_asin_tau_DL<-model_metareg_tau_spline_asin$pred$se_DL_rcs_tau
  meta1$pred_rcs_asin_tau_REML<-model_metareg_tau_spline_asin$pred$pred_REML_rcs_tau
  meta1$se_rcs_asin_tau_REML<-model_metareg_tau_spline_asin$pred$se_REML_rcs_tau
  
 
  #####################################################
  #### fractional polynomials (untransformed C-index)
  model_mfp_lme<-metamean_mfp(meta1, co=co80)
  results[nn:(nn+1),]<-model_mfp_lme$res; nn<-nn+2
  meta1$pred_fracPoly_REML<-model_mfp_lme$pred$pred_REML
  ####---------------------------------------------
  #### fractional polynomials (logit-transformed C-index)
  model_mfp_lme_logit<-metamean_mfp(meta1, cindex = "logitcindex", SEcindex = "SElogitcindex",
                                    co=co80)
  results[nn:(nn+1),]<-model_mfp_lme_logit$res; nn<-nn+2
  meta1$pred_fracPoly_logit_REML<-model_mfp_lme_logit$pred$pred_REML
  ####---------------------------------------------
  #### fractional polynomials (asin-transformed C-index)
  model_mfp_lme_asin<-metamean_mfp(meta1, cindex = "asincindex",
                                   SEcindex = "SEasincindex", co=co80)
  results[nn:(nn+1),]<-model_mfp_lme_asin$res; nn<-nn+2
  meta1$pred_fracPoly_asin_REML<-model_mfp_lme_asin$pred$pred_REML
  
  
  #####################################################
  #### exponential decay meta-regression (untransformed C-index)
  model_ecpdecay_nlme<-metamean_expdecay_lw3_R0free(meta1, co=co80)
  results[nn,]<-model_ecpdecay_nlme$res; nn<-nn+1
  meta1$pred_expdecay_REML<-model_ecpdecay_nlme$pred$pred_exponentialdecay_REML
  meta1$pred_expdecay_asymp_REML<-model_ecpdecay_nlme$pred$pred_exponentialdecay_asymp_REML
  ####---------------------------------------------
  ####  exponential decay meta-regression (logit-transformed C-index)
  model_expdecay_nlme_logit<-metamean_expdecay_lw3_R0free(meta1, cindex = "logitcindex",
                                                           invVar = "invVar_logit", co=co80)
  results[nn,]<-model_expdecay_nlme_logit$res; nn<-nn+1
  meta1$pred_expdecay_logit_REML<-model_expdecay_nlme_logit$pred$pred_exponentialdecay_REML
  meta1$pred_expdecay_logit_asymp_REML<-model_expdecay_nlme_logit$pred$pred_exponentialdecay_asymp_REML
  ####---------------------------------------------
  #### exponential decay meta-regression (asin-transformed C-index)
  model_expdecay_nlme_asin<-metamean_expdecay_lw3_R0free(meta1, cindex = "asincindex",
                                                          invVar = "invVar_asin", co=co80)
  results[nn,]<-model_expdecay_nlme_asin$res; nn<-nn+1
  meta1$pred_expdecay_asin_REML<-model_expdecay_nlme_asin$pred$pred_exponentialdecay_REML
  meta1$pred_expdecay_asin_asymp_REML<-model_expdecay_nlme_asin$pred$pred_exponentialdecay_asymp_REML
  
  
  ####################################################
  #### return results
  rownames(results)<-1:nrow(results)

  settings<-c(output$setting[ww,], nstudies_per_meta=nstudies_per_meta, ww=ww)
  
  if(any(!is.na(results$pred_80) & (results$pred_80==-9))){
    results[!is.na(results$pred_80) & (results$pred_80==-9), 5:15]<-NA
  }
  res<-list(results=results, meta1=meta1, settings=settings)
  
  return(res)
}


#### one simulation run:
res<-cind_meta(nstudies_per_meta=30, ww=1, file="simulation/Simulation_generate_meta_data_1.rda") 


#### all simulation runs
####---------------------------------------------------------------------------------------------
#### K=30 studies
r <- makeRegistry(file.dir="../scratch/MetaC/Simulation_run_metaanalysis_n30_1")
args <- expand.grid(ww=1:3000,
                    nstudies_per_meta=30, 
                    file="simulation/Simulation_generate_meta_data_1.rda")
ids <- batchMap(cind_meta, args=args, reg=r)
submitJobs(ids, reg = r, resources = list(walltime = 43000, partition="batch"))

r <- makeRegistry(file.dir="../scratch/MetaC/Simulation_run_metaanalysis_n30_2")
args <- expand.grid(ww=1:3000,
                    nstudies_per_meta=30, 
                    file="simulation/Simulation_generate_meta_data_2.rda")
ids <- batchMap(cind_meta, args=args, reg=r)
submitJobs(ids, reg = r, resources = list(walltime = 43000, partition="batch"))

r <- makeRegistry(file.dir="../scratch/MetaC/Simulation_run_metaanalysis_n30_3")
args <- expand.grid(ww=1:3000,
                    nstudies_per_meta=30, 
                    file="simulation/Simulation_generate_meta_data_3.rda")
ids <- batchMap(cind_meta, args=args, reg=r)
submitJobs(ids, reg = r, resources = list(walltime = 43000, partition="batch"))


####---------------------------------------------------------------------------------------------
#### K=50 studies
r <- makeRegistry(file.dir="../scratch/MetaC/Simulation_run_metaanalysis_n50_1")
args <- expand.grid(ww=1:3000,
                    nstudies_per_meta=50, 
                    file="simulation/Simulation_generate_meta_data_1.rda")
ids <- batchMap(cind_meta, args=args, reg=r)
submitJobs(ids, reg = r, resources = list(walltime = 43000, partition="batch"))
dim(ids)

r <- makeRegistry(file.dir="../scratch/MetaC/Simulation_run_metaanalysis_n50_2") 
args <- expand.grid(ww=1:3000,
                    nstudies_per_meta=50, 
                    file="simulation/Simulation_generate_meta_data_2.rda")
ids <- batchMap(cind_meta, args=args, reg=r)
submitJobs(ids, reg = r, resources = list(walltime = 43000, partition="batch"))

r <- makeRegistry(file.dir="../scratch/MetaC/Simulation_run_metaanalysis_n50_3")
args <- expand.grid(ww=1:3000,
                    nstudies_per_meta=50, 
                    file="simulation/Simulation_generate_meta_data_3.rda")
ids <- batchMap(cind_meta, args=args, reg=r)
submitJobs(ids, reg = r, resources = list(walltime = 43000, partition="batch"))


####---------------------------------------------------------------------------------------------
#### K=15 studies
r <- makeRegistry(file.dir="../scratch/MetaC/Simulation_run_metaanalysis_n15_1")
args <- expand.grid(ww=1:3000,
                    nstudies_per_meta=15, 
                    file="simulation/Simulation_generate_meta_data_1.rda")
ids <- batchMap(cind_meta, args=args, reg=r)
submitJobs(ids, reg = r, resources = list(walltime = 43000, partition="batch"))

r <- makeRegistry(file.dir="../scratch/MetaC/Simulation_run_metaanalysis_n15_2")
args <- expand.grid(ww=1:3000,
                    nstudies_per_meta=15, 
                    file="simulation/Simulation_generate_meta_data_2.rda")
ids <- batchMap(cind_meta, args=args, reg=r)
submitJobs(ids, reg = r, resources = list(walltime = 43000, partition="batch"))

r <- makeRegistry(file.dir="../scratch/MetaC/Simulation_run_metaanalysis_n15_3")
args <- expand.grid(ww=1:3000,
                    nstudies_per_meta=15, 
                    file="simulation/Simulation_generate_meta_data_3.rda")
ids <- batchMap(cind_meta, args=args, reg=r)
submitJobs(ids, reg = r, resources = list(walltime = 43000, partition="batch"))




####---------------------------------------------------------------------------------------------
#### when simulation runs done, aggregate results

niceResults<-function(results){
  
  settings<-data.frame(matrix(unlist(lapply(results, function(x){x$settings})), 
                              ncol=length(results[[1]]$settings), nrow=length(results), 
                              byrow=T))
  colnames(settings)<-names(results[[1]]$settings)
  
  dsdRE<-unique(settings$sdRE)
  nsm<-unique(settings$nstudies_per_meta)
  beta<-unique(settings$beta)
  tmax<-unique(settings$taumax)
  
  for(j1 in 1:length(dsdRE)){
    sdRE<-dsdRE[j1]
    
    setting<-settings[settings$taumax==tmax & settings$beta==beta & settings$sdRE==sdRE & settings$nstudies_per_meta==nsm,]
    w1<-which(settings$nstudies_per_meta==nsm & settings$sdRE==sdRE & settings$beta==beta & settings$taumax==tmax)
    
    results_list<-results_detailed<-ErrorWarningMessages<-list()
    
    for(k in 1:length(w1)){
      index<-w1[k]
      results_list[[k]]<-results[[index]]$results[,-2]
      results_list[[k]]$simrun<-k
      
      results_detailed[[k]]<-results[[index]]$meta1
    }
    
    res<-do.call("rbind", results_list)	
    
    output<-list(results=res, results_detailed=results_detailed, 
                 setting=setting)
    
    save(output, file=paste0("niceResults_Paper_Rcode/Simulation_results_sdRE",sdRE,"_beta",beta,"_tmax",tmax,"_nstudies",nsm,"_Auswertung.rda"))
    
    
  }
}


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r<-loadRegistry("../scratch/MetaC/Simulation_run_metaanalysis_n30_2")
getStatus(reg=r)
results<-reduceResultsList(reg=r)
niceResults(results)

r<-loadRegistry("../scratch/MetaC/Simulation_run_metaanalysis_n30_2")
getStatus(reg=r)
results<-reduceResultsList(reg=r)
niceResults(results)

r<-loadRegistry("../scratch/MetaC/Simulation_run_metaanalysis_n30_3")
getStatus(reg=r)
results<-reduceResultsList(reg=r)
niceResults(results)

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r<-loadRegistry("../scratch/MetaC/Simulation_run_metaanalysis_n50_1")
getStatus(reg=r)
results<-reduceResultsList(reg=r)
niceResults(results)

r<-loadRegistry("../scratch/MetaC/Simulation_run_metaanalysis_n50_2")
getStatus(reg=r)
results<-reduceResultsList(reg=r)
niceResults(results)

r<-loadRegistry("../scratch/MetaC/Simulation_run_metaanalysis_n50_3")
getStatus(reg=r)
results<-reduceResultsList(reg=r)
niceResults(results)

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r<-loadRegistry("../scratch/MetaC/Simulation_run_metaanalysis_n15_1")
getStatus(reg=r)
results<-reduceResultsList(reg=r)
niceResults(results)

r<-loadRegistry("../scratch/MetaC/Simulation_run_metaanalysis_n15_2") 
getStatus(reg=r)
results<-reduceResultsList(reg=r)
niceResults(results)

r<-loadRegistry("../scratch/MetaC/Simulation_run_metaanalysis_n15_3")
getStatus(reg=r)
results<-reduceResultsList(reg=r)
niceResults(results)















