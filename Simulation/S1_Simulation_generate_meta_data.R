library("batchtools")

### function to generate data of nstudies for meta-analysis
cind_function<-function(beta, nmax_samplesize, nstudies, seed, 
                        mu=0, sd1=0.5, gumbelLoc=0, gumbelsc=0.2, censrate=1, sc=1.1,
                        taumin=0.1, taumax=2, sdRE=0.1){
  
  library(survival)
  library(pec)
  library(survC1)
  source("S0_functions.R") 
  
  ####--------------------------------------------------------------------------
  #### prepare empty data frames to store results
  cindices_bootstrap<-list()
  cind<-matrix(NA,nrow=nstudies, ncol=12)
  colnames(cind)<-c("id", "n", "nevents", "tau", "study_beta",
                    "cindex", "SEcindex", 
                    "logitcindex", "SElogitcindex", 
                    "asincindex", "SEasincindex",
                    "true_CensoringProp")
  
  #### set seed
  set.seed(seed)
  ####--------------------------------------------------------------------------
  ### draw study follow-up times 
  taus<-qgamma(seq(0.01, 0.95, 0.001), shape=1.5, rate=1)
  taus<-taus[taus>taumin & taus<taumax] #### maximal follow-up time
  
  #### ensure, that at least one study's FU longer than cut-off
  studytau<-0
  while(max(studytau)<0.8*taumax){
    ## sample administrative censoring times per study
    studytau<-taus[sample(1:length(taus), nstudies, replace=T)] 
  }
  
  ### draw study-specific random effect:
  reff<-rnorm(nstudies, mean=0, sd = sdRE) 
  
  #### for each study b in nstudies:    
  for(b in 1:nstudies){
    print(b)
    set.seed(1000*seed+b)
    
    logc<-c2<-selogc<-cilow<-ciup<-NA
    count<-1
    nevents<-0
    
    while(nevents<10){ #### ensure each study has at least 10 events
      
      if(count>1){
        ## if FU too short (for 10 events) draw another tau
        studytau[b] <-taus[sample(1:length(taus), 1)] 
      }
      
      ## counter to avoid never-ending loops
      count2<-1
      while(nevents<10 & count2<100){ 
        #### draw sample size
        ntemp<-sample(seq(100,nmax_samplesize,10), 1)
        #### Weibull distributed survival times
        dat<-getWeibullData(n=ntemp, gumbelsc = gumbelsc, gumbelLoc = gumbelLoc, mu =mu,
                            seed = 1000*seed+b+count,sd1 = sd1, sc = sc , beta =beta, 
                            censrate = censrate,silent = T)$data
        #### censoring status and censoring time
        dat$status_new <- ((dat$Ttrue <= dat$Ctrue) & dat$Ttrue<studytau[b]) + 0 # 
        dat$Ttilde_new <- pmin(dat$Ttrue, dat$Ctrue, studytau[b])
        #### number of events
        nevents<-sum(dat$status_new==1)
        count<-count+1
        count2<-count2+1
      }
    }
    
    
    ####-------------------------------    
    ## add random effect:
    tc <- coxph(Surv(Ttilde_new,status_new)~X1,data=dat, x=T)
    cindex_pec <- cindex(tc, formula = Surv(Ttilde_new,status_new)~X1, data=dat)$AppCindex$coxph
    cindex_pec<-cindex_pec+reff[b]

    ####  bootstrap variance of C-index per study
   ctemp<-vector()
   for(kk in 1:1000){
    nevents<-0
    while(nevents<5){
      temp<-dat[sample(1:nrow(dat), nrow(dat), replace = T),]  
      nevents<-sum(temp$status_new==1)
    }
    tctemp<- try(coxph(Surv(Ttilde_new,status_new)~X1,data=temp, x=T))
    try(ctemp[kk]<-cindex(tctemp, formula = Surv(Ttilde_new,status_new)~X1, data=temp)$AppCindex$coxph)
    try(ctemp[kk]<-ctemp[kk]+reff[b])
            
  }
          
  #### SE of untransformed C-Index
  try(SEcindex<-sd(ctemp, na.rm=T))
  #### transformed C-Indices
  try(logitc<-log(cindex_pec/(1-cindex_pec)))
  try(asinc<-asin(sqrt(cindex_pec)))
          
  ##### SE of transformed C-Indices
  try(selogitc<-sd(log(ctemp/(1-ctemp)), na.rm=T))
  try(seasinc<-sd(asin(sqrt(ctemp)), na.rm=T))
          
  try(cind[b,]<-c(b, ntemp, sum(dat$status_new==1) ,studytau[b], beta,
                  cindex_pec,  SEcindex,  logitc, selogitc,
                  asinc, seasinc, 1-sum(dat$status_new)/nrow(dat)))
  }
        
  cind<-data.frame(cind)
        
  #### store settings
  settings<-list(settings=c("beta"=beta, "mu"=mu, "sd"=sd1, "gumbelLoc"=gumbelLoc, 
                            "gumbelsc"=gumbelsc, "censrate"=censrate, "sc"=sc, 
                            "nstudies"=nstudies, "nmax_samplesize"=nmax_samplesize,
                            "sdRE"=sdRE,"taumin"=taumin, "taumax"=taumax),
                 taus=taus)
  
    
  results<-list("studies"=cind, "settings"=settings) 
        
  return(results)
}
    

### try one simulation run:
ex<-cind_function(beta=0.5, nmax_samplesize=300, nstudies=1, seed=1, 
  mu=0, sd1=0.5, gumbelLoc=0, gumbelsc=0.2, censrate=0.5, sc=1.1,
  taumin=0.1, sdRE=0)
    
    
################################################################################
########### multiple simulation runs for different heterogeneities and follow-up times
r <- makeRegistry(file.dir="../scratch/MetaC/Simulation_generate_meta_data_1")
args <- expand.grid(seed=1001:2000,
                    nmax_samplesize=1000,
                    nstudies=50,
                    beta=0.5,
                    censrate=0.5,
                    sdRE=c(0, 0.01,0.03),
                    taumax=0.7) 
ids <- batchMap(cind_function, args=args, reg=r)
submitJobs(ids, reg = r, resources = list(walltime = 43000, partition="batch"))
  
  
########### Setting
r <- makeRegistry(file.dir="../scratch/MetaC/Simulation_generate_meta_data_2")
args <- expand.grid(seed=2001:3000,
                    nmax_samplesize=1000,
                    nstudies=50,
                    beta=0.5,
                    censrate=0.5,
                    sdRE=c(0, 0.01,0.03),
                    taumax=2) 
ids <- batchMap(cind_function, args=args, reg=r)
submitJobs(ids, reg = r, resources = list(walltime = 43000, partition="batch"))
    
    
########### Setting 
r <- makeRegistry(file.dir="../scratch/MetaC/Simulation_generate_meta_data_3")
args <- expand.grid(seed=3001:4000,
                    nmax_samplesize=1000,
                    nstudies=50,
                    beta=0.5,
                    censrate=0.5,
                    sdRE=c(0, 0.01,0.03),
                    taumax=0.9) 
ids <- batchMap(cind_function, args=args, reg=r)
submitJobs(ids, reg = r, resources = list(walltime = 43000, partition="batch"))    
    


################################################################################
#### when simulation runs done; aggregate results:
niceResults<-function(results){
  
  studies<-lapply(results, function(x){x$studies})
  setting<-matrix(unlist(lapply(results, function(x){x$settings$settings})), ncol=length(results[[1]]$settings$settings),nrow=length(results), byrow=T)
  colnames(setting)<-names(results[[1]]$settings$settings)
  
  return(list(studies=studies, setting=setting))
}

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

r<-loadRegistry("../scratch/MetaC/Simulation_generate_meta_data_1")
getStatus(reg=r)
results<-reduceResultsList(reg=r)
output<-niceResults(results)
save(output, file="data/Simulation_generate_meta_data_1.rda")

r<-loadRegistry("../scratch/MetaC/Simulation_generate_meta_data_2")
getStatus(reg=r)
results<-reduceResultsList(reg=r)
output<-niceResults(results)
save(output, file="data/Simulation_generate_meta_data_2.rda")

r<-loadRegistry("../scratch/MetaC/Simulation_generate_meta_data_3")
getStatus(reg=r)
results<-reduceResultsList(reg=r)
output<-niceResults(results)
save(output, file="data/Simulation_generate_meta_data_3.rda")



    
