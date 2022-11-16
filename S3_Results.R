## function to generate result figures and tables

library(ggplot2)
library(gridExtra)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(purrr)


#### different settings
sdRE<-c(0, 0.01, 0.03) 
taumax<-c(0.7,0.9, 2)


for(nstudies in c(15, 30, 50)){
  #### list of files with nstudies
  loadnames<-paste0("simulation/Simulation_results_sdRE",
                    rep(sdRE, each=length(taumax)), "_beta0.5_tmax", 
                    taumax, "_nstudies",nstudies,"_Auswertung.rda")
  #### define prefix for storing figures and tables
  savename<-paste0("results/MetaC_nstudies", nstudies)

  results_list<-settings<-list()

  for(nn in seq_along(loadnames)){
    #### load files
    load(loadnames[nn])
    #### remove (accidental) spaces from names
    output$results$name<-gsub(" ", "", output$results$name)
    
    #### store results in list
    results_list[[nn]]<-output$results

    #### store all settings
    settings[[nn]]<-output$setting[!duplicated(output$setting[,-14]),]
    results_list[[nn]]$beta<-settings[[nn]]$beta
    results_list[[nn]]$mu<-settings[[nn]]$mu
    results_list[[nn]]$sd<-settings[[nn]]$sd
    results_list[[nn]]$gumbelLoc<-settings[[nn]]$gumbelLoc
    results_list[[nn]]$gumbelsc<-settings[[nn]]$gumbelsc
    results_list[[nn]]$censrate<-settings[[nn]]$censrate
    results_list[[nn]]$sc<-settings[[nn]]$sc
    results_list[[nn]]$nstudies<-settings[[nn]]$nstudies
    results_list[[nn]]$nmax_samplesize<-settings[[nn]]$nmax_samplesize
    results_list[[nn]]$sdRE<-settings[[nn]]$sdRE
    results_list[[nn]]$taumin<-settings[[nn]]$taumin
    results_list[[nn]]$taumax<-settings[[nn]]$taumax
    results_list[[nn]]$nstudies_per_meta<-settings[[nn]]$nstudies_per_meta
    rm(output)
    }
    
    ### generate one large data-frame
    results <- map_df(results_list, ~as.data.frame(.x))
    results$var_method_LW<-results$var_method[1:39]
    results$id<-paste0(results$name, "_", results$var_method_LW)
    rownames(results)<-1:nrow(results)
    
    ####------------------------------------------------------------------------
    #### store settings as data.frame
    settings<-results[,c(19:31)]
    #### extract settings with nstudies
    whichSetting<-which(settings$nstudies_per_meta==nstudies)
    settings<-settings[whichSetting,!(colnames(settings) %in% "ww")]

    
    #### vector of nices names for all methods
    nicenames<-c('raw_REML'= 'MA (id)' ,
                 'raw_last50perc_REML'= 'MA (id, last 50 %)',
                 'raw_last30perc_REML'= 'MA (id, last 30 %)',
                 'logit_REML'= 'MA (logit)',
                 'logit_last50perc_REML'= 'MA (logit, last 50 %)', 
                 'logit_last30perc_REML'= 'MA (logit, last 30 %)', 
                 'asin_REML'= 'MA (asin)', 
                 'asin_last50perc_REML'= 'MA (asin, last 50 %)', 
                 'asin_last30perc_REML'= 'MA (asin, last 30 %)',
                 'raw_tau_REML' ='linear (id)', 
                 'logit_tau_REML'='linear (logit)', 
                 'asin_tau_REML'='linear (asin)', 
                 'rmaspline_cindex_REML'='RCS (id)',
                 'rmaspline_logitcindex_REML'='RCS (logit)',
                 'rmaspline_asincindex_REML'='RCS (asin)',
                 'fracpoly_cindex_REML'='FP2 (id)', 
                 'fracpoly_logitcindex_REML'='FP2 (logit)',
                 'fracpoly_asincindex_REML'='FP2 (asin)',
                 'raw_DL'= 'raw (DL)', 
                 'raw_last50perc_DL'= 'raw (last 50 %, DL)',
                 'raw_last30perc_DL'= 'raw (last 30 %, DL)',
                 'logit_DL'= 'logit (DL)', 
                 'logit_last50perc_DL'= 'logit (last 50 %, DL)',
                 'logit_last30perc_DL'= 'logit (last 30 %, DL)',              
                 'asin_DL'= 'asin (DL)',
                 'asin_last50perc_DL'= 'asin (last 50 %, DL)',
                 'asin_last30perc_DL'= 'asin (last 30 %, DL)',
                 'raw_tau_DL' ='linear (raw, DL)',
                 'logit_tau_DL' ='linear (logit, DL)', 
                 'asin_tau_DL' ='linear (asin, DL)',
                 'rmaspline_cindex_DL'='RCS (DL)',
                 'rmaspline_logitcindex_DL'='RCS (logit, DL)',
                 'rmaspline_asincindex_DL'='RCS (asin, DL)',
                 'fracpoly_cindex_DL'='frac. poly. (DL)', 
                 'fracpoly_logitcindex_DL'='frac. poly. (logit, DL)',
                 'fracpoly_asincindex_DL'='frac. poly. (asin, DL)',
                 'expdecay_cindex_REML'='exp. decay (REML)',
                 'expdecay_logitcindex_REML'='exp. decay (logit, REML)', 
                 'expdecay_asincindex_REML'='exp. decay (asin, REML)')
    
    #### match names to get nice names
    results$nicenames<-nicenames[match(results$id, names(nicenames))]
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### figures preparation
    tempResults<-results
    #### prep data for results tables and figures
    tempResults$method<-tempResults$name
    
    namefuersdRE<-"~(sigma[a]=="
    namefuerTaumax<-"~(tau[max]=="
    tempResults$label_sdRE<- factor(NA,
                         levels=c(paste0("no~heterogeneity", namefuersdRE,"0)"),
                                  paste0("moderate~heterogeneity", namefuersdRE,"0.01)"),
                                  paste0("high~heterogeneity", namefuersdRE,"0.03)")))
    tempResults$label_sdRE[tempResults$sdRE==0]<-paste0("no~heterogeneity", namefuersdRE,tempResults$sdRE[tempResults$sdRE==0], ")")
    tempResults$label_sdRE[tempResults$sdRE==0.01]<-paste0("moderate~heterogeneity", namefuersdRE,tempResults$sdRE[tempResults$sdRE==0.01], ")")
    tempResults$label_sdRE[tempResults$sdRE==0.03]<-paste0("high~heterogeneity", namefuersdRE,tempResults$sdRE[tempResults$sdRE==0.03], ")")
    
    tempResults$label_taumax<- factor(NA,
                           levels=c(paste0("short~FU", namefuerTaumax,"0.7)"),
                                    paste0("moderate~FU", namefuerTaumax,"0.9)"),
                                    paste0("long~FU", namefuerTaumax,"2)")))
    tempResults$label_taumax[tempResults$taumax==0.7]<-paste0("short~FU", namefuerTaumax,tempResults$taumax[tempResults$taumax==0.7], ")")
    tempResults$label_taumax[tempResults$taumax==0.9]<-paste0("moderate~FU", namefuerTaumax,tempResults$taumax[tempResults$taumax==0.9], ")")
    tempResults$label_taumax[tempResults$taumax==2]<-paste0("long~FU", namefuerTaumax,tempResults$taumax[tempResults$taumax==2], ")")
    
    tempResults$label_sdRE<-factor( tempResults$label_sdRE)
    tempResults$label_taumax<-factor( tempResults$label_taumax)
    
    tempResults$name<-paste0("sdRE", tempResults$sdRE, "_beta", 
                             tempResults$beta, "_tmax", tempResults$taumax)
    
  
    tempErrorWarning<-tempResults
    ####------------------------------------------------------------------------
    #### results tables
    tempResults<-tempResults[!is.na(tempResults$pred_80),]
    #### prediction for 80% of maximal follow-up time
    agresults_mean<-aggregate(tempResults$pred_80,
                              by=list(tempResults$id, tempResults$name), mean)
    tabresults1<-reshape2::dcast(agresults_mean, 
                                 Group.1 ~  Group.2, value.var = c("x"))
    agresults_sd<-aggregate(tempResults$pred_80, 
                            by=list(tempResults$id, tempResults$name), sd)
    tabresults2<-reshape2::dcast(agresults_sd, 
                                 Group.1 ~  Group.2, value.var = c("x"))
    #### area between predicted curve and true curve
    AbCresults_mean<-aggregate(tempResults$AbC,
                               by=list(tempResults$id, tempResults$name), mean)
    tabresults3<-reshape2::dcast(AbCresults_mean, 
                                 Group.1 ~  Group.2, value.var = c("x"))
    AbCresults_sd<-aggregate(tempResults$AbC, 
                             by=list(tempResults$id, tempResults$name), sd)
    tabresults4<-reshape2::dcast(AbCresults_sd, Group.1 ~  Group.2, 
                                 value.var = c("x"))

    #### create nice output tables
    settings<-unique(tempResults$name)
    tablist_pred<-tablist_abc<-list()
    for(ss in seq_along(settings)){
      sett<-settings[ss]
      #### prediction for 80% of maximal follow-up time
      m1<-merge(tabresults1[,c("Group.1", sett)], 
                tabresults2[,c("Group.1", sett)], by="Group.1", 
                suffixes = c("_mean", "_sd"))
      m1$nicenames<-nicenames[match(m1$Group.1, names(nicenames))]
      
      m1<-m1[,c(4,2,3)]
      colnames(m1)<-c("method", "mean", "sd")
      
      m1$method<-factor(m1$method, levels=nicenames)
      m1<-m1[match(levels(m1$method),m1$method),]
      tablist_pred[[ss]]<-m1
      rm(m1)
      
      
      #### area between predicted curve and true curve
      m1<-merge(tabresults3[,c("Group.1", sett)], 
                tabresults4[,c("Group.1", sett)], by="Group.1", 
                suffixes = c("_mean", "_sd"))
      m1$nicenames<-nicenames[match(m1$Group.1, names(nicenames))]
      
      m1<-m1[,c(4,2,3)]
      colnames(m1)<-c("method", "mean", "sd")
      
      m1$method<-factor(m1$method, levels=nicenames)
      m1<-m1[match(levels(m1$method),m1$method),]
      tablist_abc[[ss]]<-m1
      rm(m1)
      
    }
    
    names(tablist_pred)<-names(tablist_abc)<-settings
    
    save(tablist_pred,  tablist_abc,
         file=paste0(savename,"_results_estimate_and_area_Table.rda"))
    
    ####------------------------------------------------------------------------
    #### FIGURES
    ## asymptotic C-index
    asympC<-0.7446737   
    # C-index at 80perc of max follow-up times
    C80<-rep(c(0.7774458, 0.7665837, 0.7448051),3)
    
    dummy2<-expand.grid(label_taumax=c('short~FU~(tau[max]==0.7)', 'moderate~FU~(tau[max]==0.9)', 'long~FU~(tau[max]==2)'),
                        label_sdRE=c('no~heterogeneity~(sigma[a]==0)', 'moderate~heterogeneity~(sigma[a]==0.01)', 'high~heterogeneity~(sigma[a]==0.03)'))
    
    dummy2$C80<-C80
    
    
    #### methods with different transformations
    regression_REML_id<-c('MA (id)',  'MA (id, last 50 %)', 
                          'MA (id, last 30 %)',  'linear (id)', 
                          'RCS (id)',  
                          'FP2 (id)')
    regression_REML_logit<-c('MA (logit)',  'MA (logit, last 50 %)', 
                             'MA (logit, last 30 %)',  'linear (logit)', 
                             'RCS (logit)',  
                             'FP2 (logit)')
    regression_REML_asin<-c('MA (asin)',  'MA (asin, last 50 %)', 
                            'MA (asin, last 30 %)',  'linear (asin)', 
                            'RCS (asin)',  
                            'FP2 (asin)')
    
    
    ####------------------------------------------------------------------------
    ##### transformation-based figures:
    #### color palette
    cbPalette <- c("#7bbcb0","#94d1c6","#b3e8df",
                   "#8da8cc", "#bc97bd", 
                   "#db9f69","#f57d7d")
    
    
    ####--------------------------------------
    #### untransformed meta-regression models
    temp2<-tempResults[tempResults$nicenames %in% regression_REML_id,]
    temp2REML<-temp2
    temp2REML$nicenames<-factor(temp2REML$nicenames, levels=regression_REML_id)
    
    
    ggplot(temp2REML, 
           aes(x=nicenames, y=pred_80, group=nicenames,
               fill=factor(nicenames))) + 
      geom_boxplot(lwd=0.2, outlier.size=0.2)+ 
      geom_hline(aes(yintercept=asympC), col="gray60")+
      geom_hline(data=dummy2, aes(yintercept=C80), col="red")+
      facet_grid(label_taumax~label_sdRE, scales="free_y", 
                 labeller = labeller(.cols=label_parsed,
                                     .rows=label_parsed,
                                     .multi_line=FALSE))+
      theme_bw()+ theme(axis.text.x = element_text(angle = 40, hjust = 1, 
                                                   color=cbPalette, size=12),
                        axis.title.y=element_text(size=14),
                        axis.title.x=element_blank(),
                        panel.grid.major.x = element_blank() ,
                        panel.grid.major.y = element_line( size=.1, 
                                                           color="gray90")) +
      guides(fill=FALSE)+ ylab("Pooled C-Index")+
      scale_fill_manual(values=cbPalette)+ylim(c(0.7,0.8))
    
    ggsave(paste0(savename,"boxplots_metaregression_IDmethods.pdf"), 
           width = 20, height=20, units = "cm")
    
    
    ####--------------------------------------
    ####logit transformed meta-regression models
    temp2<-tempResults[tempResults$nicenames %in% regression_REML_logit,]
    temp2REML<-temp2
    temp2REML$nicenames<-factor(temp2REML$nicenames, 
                                levels=regression_REML_logit)
    
    ggplot(temp2REML, 
           aes(x=nicenames, y=pred_80, group=nicenames, 
               fill=factor(nicenames))) + 
      geom_boxplot(lwd=0.2, outlier.size=0.2)+ 
      geom_hline(aes(yintercept=asympC), col="gray60")+
      geom_hline(data=dummy2, aes(yintercept=C80), col="red")+
      facet_grid(label_taumax~label_sdRE, scales="free_y", 
                 labeller = labeller(.cols=label_parsed,
                                     .rows=label_parsed,
                                     .multi_line=FALSE))+ 
      theme_bw()+ theme(axis.text.x = element_text(angle = 40, hjust = 1, 
                                                   color=cbPalette, size=12),
                        axis.title.y=element_text(size=14),
                        axis.title.x=element_blank(),
                        panel.grid.major.x = element_blank() ,
                        panel.grid.major.y = element_line(size=.1, 
                                                          color="gray90" ) ) +
      guides(fill=FALSE)+ ylab("Pooled C-Index")+
      scale_fill_manual(values=cbPalette)+ylim(c(0.7,0.8))
    
    ggsave(paste0(savename,"boxplots_metaregression_logitmethods.pdf"), 
           width = 20, height=20, units = "cm")
    
    ####--------------------------------------
    #### asin transformed meta-regression models
    temp2<-tempResults[tempResults$nicenames %in% regression_REML_asin,]
    temp2REML<-temp2
    temp2REML$nicenames<-factor(temp2REML$nicenames, 
                                levels=regression_REML_asin)
    
    ggplot(temp2REML, 
           aes(x=nicenames, y=pred_80, group=nicenames, 
               fill=factor(nicenames))) + 
      geom_boxplot(lwd=0.2, outlier.size=0.2)+ 
      geom_hline(aes(yintercept=asympC), col="gray60")+
      geom_hline(data=dummy2, aes(yintercept=C80), col="red")+
      facet_grid(label_taumax~label_sdRE, scales="free_y", 
                 labeller = labeller(.cols=label_parsed,
                                     .rows=label_parsed,
                                     .multi_line=FALSE))+ 
      theme_bw()+ theme(axis.text.x = element_text(angle = 40, hjust = 1, 
                                                   color=cbPalette, size=12),
                        axis.title.y=element_text(size=14),
                        axis.title.x=element_blank(),
                        panel.grid.major.x = element_blank() ,
                        panel.grid.major.y = element_line( size=.1, 
                                                           color="gray90")) +
      guides(fill=FALSE)+ ylab("Pooled C-Index")+
      scale_fill_manual(values=cbPalette)+ylim(c(0.7,0.8))
    
    ggsave(paste0(savename,"boxplots_metaregression_asinmethods.pdf"), 
           width = 20, height=20, units = "cm")
    
    
    ############################################################################
    ####------------------------------------------------------------------------
    #### FAILURES
    tempErrorWarning<-tempErrorWarning[grepl("REML", 
                                             tempErrorWarning$id),] #### only REML methods
    tempErrorWarning$errwarn<-((tempErrorWarning$anyError + 
                                  tempErrorWarning$anyWarning)>0)*1
    tempErrorWarning$errwarn[is.na(tempErrorWarning$pred_80)]<-1
    temp3<-aggregate(tempErrorWarning$errwarn, 
                     by=list("id"=tempErrorWarning$id, 
                             label_sdRE=tempErrorWarning$label_sdRE, 
                             label_taumax=tempErrorWarning$label_taumax), mean)
    temp3$nErrorsWarnings<-aggregate(tempErrorWarning$errwarn,
                                     by=list("id"=tempErrorWarning$id, 
                                             label_sdRE=tempErrorWarning$label_sdRE,
                                             label_taumax=tempErrorWarning$label_taumax), sum)$x
    temp3$x<-temp3$x*100
    temp3$nicenames<-nicenames[match(temp3$id, names(nicenames))]
    temp3$nicenames<-factor(temp3$nicenames, 
                            levels=nicenames[nicenames%in%temp3$nicenames])
    
    #### FAILURES - Figures
    ggplot(temp3, aes(x=nicenames, y=x ))+
      geom_bar(stat="identity", position="stack")+
      # scale_fill_manual(values = c("darkgoldenrod1", "brown2"))+ 
      facet_grid(label_taumax~label_sdRE, 
                 labeller = labeller(.cols=label_parsed,
                                     .rows=label_parsed,
                                     .multi_line=FALSE))+ ylim(0,100)+
      theme_bw()+ theme(axis.text.x = element_text(angle = 40, hjust = 1, size=6),
                        axis.title.x=element_blank(),
                        panel.grid.major.x = element_blank() ,
                        panel.grid.major.y = element_line( size=.1, color="gray90" ) ) +
      guides(fill=FALSE)+ ylab("rel. frequency of errors (red) and warnings (yellow)") 
    ggsave(paste0(savename,"_barplot_nErrorsWarnings.pdf"), width = 20, height=20, units = "cm")
    
    ####------------------------------------------------------------------------------------
    #### FAILURES - Table
    temp3<-temp3[,c( "nicenames","nErrorsWarnings","label_sdRE","label_taumax",  "x")]
    temp3$setting<-paste0("sdRE",gsub(".*= ","",temp3$label_sdRE),"_beta0.5_tmax",
                          gsub(".*= ","",temp3$label_taumax))
    
    setDT(temp3)
    settings<-unique(temp3$setting)
    
    errwarnlist<-list()
    nruns<-1000
    for(ss in seq_along(settings)){
      temp5<-temp3[temp3$setting==settings[ss],]
      temp5<-temp5[,c('nicenames','x')]
      
      errwarnlist[[ss]]<-temp5
    }
    
    names(errwarnlist)<-settings
    temp<-errwarnlist
    for(i in 1:9){
      colnames(temp[[i]])[2]<-"propErrWarn"
    }
    

    names(temp)[c(1,4,7,2,5,8,3,6,9)]
    papertable_failure<-cbind(temp[[1]], temp[[4]][,2], temp[[7]][,2],
                              temp[[2]][,2], temp[[5]][,2], temp[[8]][,2],
                              temp[[3]][,2], temp[[6]][,2], temp[[9]][,2])
    colnames(papertable_failure)<-
      c("method", 
        "short_FU_sdRE0", "moderate_FU_sdRE0", "long_FU_sdRE0",
        "short_FU_sdRE0.01", "moderate_FU_sdRE0.01", "long_FU_sdRE0.01", 
        "short_FU_sdRE0.03", "moderate_FU_sdRE0.03", "long_FU_sdRE0.03")
    

    papertable_failure$nicenames<-factor(papertable_failure$method, 
                               levels=nicenames[nicenames %in% 
                                                  papertable_failure$method])
    
    save(papertable_failure,
         file=paste0(savename,"_Table_nErrorsWarnings.rda"))
    ####------------------------------------------------------------------------
    #### COVERAGE
    cbPalette <- c("#7bbcb0","#94d1c6","#b3e8df",
                   "#8da8cc", "#bc97bd", 
                   "#db9f69","#f57d7d")
    
    
    #### Table of coverage rate
    temp3<-aggregate(tempErrorWarning$inCI, 
                     by=list("id"=tempErrorWarning$id, 
                             label_sdRE=tempErrorWarning$label_sdRE, 
                             label_taumax=tempErrorWarning$label_taumax), mean,
                     na.rm=T)
    temp3<-temp3[!is.na(temp3$x),]
    temp3$x<-temp3$x*100
    temp3$nicenames<-nicenames[match(temp3$id, names(nicenames))]
    temp3$nicenames<-
      factor(temp3$nicenames, levels=nicenames[nicenames%in%temp3$nicenames])
    temp3$se<-qt(0.975, df = 999)* sqrt((temp3$x/100)*(1-temp3$x/100)/1000)*100
    
    ####---- coverage of untransformed C-index
    temp3id<-temp3[!grepl("logit", temp3$nicenames) & 
                     !grepl("asin", temp3$nicenames) & 
                     !grepl("log", temp3$nicenames),]
    temp3id<-temp3id[!grepl("MicMen|decay", temp3id$nicenames),]
    temp3id$nicenames<-factor(temp3id$nicenames)
    
    ggplot(temp3id, aes(x=nicenames, y=x))+
      geom_point()+
      geom_hline(yintercept=95, col="red")+
      facet_grid(label_taumax~label_sdRE, 
                 labeller = labeller(.cols=label_parsed,
                                     .rows=label_parsed,
                                     .multi_line=FALSE))+
      geom_errorbar(aes(ymin=x-se, ymax=x+se), width=.7,
                    position=position_dodge(0.05))+
      theme_bw()+ theme(axis.text.x = element_text(angle = 40, 
                                                   hjust = 1, size=10),
                        axis.title.y = element_text(size=12),
                        axis.title.x=element_blank(),
                        panel.grid.major.x = element_blank() ,
                        panel.grid.major.y = element_line( size=.1, 
                                                           color="gray90" ) ) +
      guides(fill=FALSE)+ ylab("Coverage [%]")

    ggsave(paste0( savename,"_Whiskerplot_CoverageRate_id.pdf"),
           width = 20, height=20, units = "cm")
    
    ####---- coverage of logit-transformed C-index
    temp3logit<-temp3[grepl("logit", temp3$nicenames),]
    temp3logit<-temp3logit[!grepl("MicMen|decay", temp3logit$nicenames),]
    temp3logit$nicenames<-factor(temp3logit$nicenames)
    
    ggplot(temp3logit, aes(x=nicenames, y=x))+
      geom_point()+
      geom_hline(yintercept=95, col="red")+
      facet_grid(label_taumax~label_sdRE, 
                 labeller = labeller(.cols=label_parsed,
                                     .rows=label_parsed,
                                     .multi_line=FALSE))+
      geom_errorbar(aes(ymin=x-se, ymax=x+se), width=.7,
                    position=position_dodge(0.05))+
      theme_bw()+ theme(axis.text.x = element_text(angle = 40, 
                                                   hjust = 1, size=10),
                        axis.title.y = element_text(size=12),
                        axis.title.x=element_blank(),
                        panel.grid.major.x = element_blank() ,
                        panel.grid.major.y = element_line( size=.1, 
                                                           color="gray90" ) ) +
      guides(fill=FALSE)+ ylab("Coverage [%]")
    ggsave(paste0(savename,"_Whiskerplot_CoverageRate_logit.pdf"), 
           width = 20, height=20, units = "cm")
    
    ####---- coverage of asin-transformed C-index
    temp3asin<-temp3[grepl("asin", temp3$nicenames),]
    temp3asin<-temp3asin[!grepl("MicMen|decay", temp3asin$nicenames),]
    temp3asin$nicenames<-factor(temp3asin$nicenames)
    
    ggplot(temp3asin, aes(x=nicenames, y=x))+
      geom_point()+
      geom_hline(yintercept=95, col="red")+
      facet_grid(label_taumax~label_sdRE, 
                 labeller = labeller(.cols=label_parsed,
                                     .rows=label_parsed,
                                     .multi_line=FALSE))+
      geom_errorbar(aes(ymin=x-se, ymax=x+se), width=.7,
                    position=position_dodge(0.05))+
      theme_bw()+ theme(axis.text.x = element_text(angle = 40, 
                                                   hjust = 1, size=10),
                        axis.title.y = element_text(size=12),
                        axis.title.x=element_blank(),
                        panel.grid.major.x = element_blank() ,
                        panel.grid.major.y = element_line( size=.1, 
                                                           color="gray90" ) ) +
      guides(fill=FALSE)+ ylab("Coverage [%]")
    ggsave(paste0(savename,"_Whiskerplot_CoverageRate_asin.pdf"),
           width = 20, height=20, units = "cm")
    

    
  }
