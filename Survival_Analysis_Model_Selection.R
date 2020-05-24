
####################################################################
###                Survial Analysis Model Selection              ###
####################################################################

library(survivalROC)
library(survival)
library(survminer)
library(survcomp)
library(ggplot2)
library(dplyr)

library(meta)
# calculate the overall HR (AUC, and C.index) for each signature, 
# based on its HR (AUC, and C.index) value and standard errors derived from individual validation sets


# TCGA-PRAD     == Cell 2015                == GDC         == RNAseq: CPM              == Time to BCR
# CPC-GENE      == Nature_2017              == GSE107299   == HuGene 2.0 ST & HTA 2.0  == Time to BCR
# MSKCC         ==_Cancer_Cell_2010         == GSE21034    == HuEx 1.0 ST              == Time to BCR
# DKFZ          == Cancer Cell 2018         == cBioPortal  == RNAseq: RPKM             == Time to BCR
# GSE54460      == Cancer Research 2014     == GSE54460    == RNAseq: CPM              == Time to BCR
# CamCap        == Nature Genetics 2016     == GSE70768    == Illumina HumanHT-12 V4.0 == Time to BCR
# Stockholm     == Nature Genetics 2016     == GSE70769    == Illumina HumanHT-12 V4.0 == Time to BCR
# DESNT         == Eur Urol Focus 2017      == GSE94767    == HuEx 1.0 ST              == Time to BCR
# E-MTAB-6128   == Annuals of Oncology 2018 == E-MTAB-6128 == HuGene 2.0 ST            == Time to BCR
# Belfast       == Journal of Clinical Oncology 2017  == GSE116918 == ADXPCv1a520642   == Time to BCR/Metastasis


#################################################################

###### Parameters

datasets <- c('TCGA_PRAD','GSE107299','GSE21034','DKFZ2018','GSE54460','GSE70768','GSE70769',
              'GSE94767','E_MTAB_6128','GSE116918_BCR') # , 'GSE116918_Metastasis'

# Published signatures
signatures <- c('Agell','Bibikova','Bismar','Decipher','Ding','Glinsky','Irshad',
                'Jennifer','Jia','Kamoun','Li','Long','Luca','Mo','Nakagawa','Olmos',
                'Oncotype','Penney','Planche','Prolaris','Ramaswamy','Ramos_Montoya',
                'Ross_Adams','Ross_Robert','Sharma','Talantov','Varambally','Wu','Yang',
                'Yu')
                
# SuperPC
models <- c('CoxPH', 'CoxNetAlpha0', 'CoxNet', 'SuperPC0', 'SuperPC0.1', 'SuperPC') # SuperPC 0.3

# plsRcox
models <- c('CoxPH', 'CoxNetAlpha0', 'CoxNet', 'plsRcox1','plsRcox2','plsRcox3') # plsRcox2

# RandomForest
models <- c('CoxPH', 'CoxNetAlpha0', 'CoxNet', 'randomForestSRC10', 'randomForestSRC20',
            'randomForestSRC50', 'randomForestSRC100') # randomForestSRC100

i <- 0
statsList <- list()

for (dataset in datasets) {
  
  message(dataset)
  
  for (model in models) {
    
    for (signature.name in signatures) {
      
      message(signature.name)
      
      i <- i + 1
      stats <- c()
      
      fl.name <- paste0('report/Survival/CV10Scale/Signature/Survival_',
                        model, '_', signature.name, '_', dataset, '.csv')
      
      # fl.name <- paste0('report/SurvivalC/CV10Scale/Signature/Survival_',
      #                  model, '_', signature.name, '_', dataset, '.csv')
      
      if (! file.exists(fl.name)) {
        print (paste('NA:', fl.name))
        stats <- c(stats, c(dataset, model, signature.name, rep(NA,10)))
        statsList[[i]] <- stats
        next
      }
      
      res <- read.csv(file=fl.name, stringsAsFactors = F, row.names = 1)
      
      if (sum(res$risk.score==0)>=length(res$risk.score)) {
        print (paste('Zero:', fl.name))
        stats <- c(stats, c(dataset, model, signature.name, rep(NA,10)))
        statsList[[i]] <- stats
        next
      }
      
      
      ###
      
      risk.score <- res$risk.score
      risk.group <- res$risk.group
      bcr.time <- res$bcr.time
      bcr.status <- res$bcr.status
      
      ############ C index
      
      c <- concordance.index(x=risk.score, 
                             surv.time=bcr.time, 
                             surv.event=bcr.status, 
                             #cl=risk.group,
                             method="noether")
      c.index <- c$c.index
      
      stats <- c(stats, c(dataset, model, signature.name, c.index))
      
      
      ############ Time-dependent ROC
      
      cutoff <- 5*12
      nobs <- length(risk.score)
      roc <- survivalROC(Stime=bcr.time,
                         status=bcr.status,
                         marker = risk.score,
                         method = 'NNE',
                         predict.time = cutoff,
                         span = 0.25*nobs^(-0.20))
      
      td.auc <- roc$AUC
      
      stats <- c(stats, td.auc)
      
      ############ CoxPH
      
      coxtest <- coxph(Surv(bcr.time, bcr.status) ~ risk.score)
      summcph <- summary(coxtest)
      
      coeffs <- c(summcph$coefficients[,2], summcph$conf.int[,3:4], 
                  summcph$coefficients[,5])
      
      stats <- c(stats, coeffs)
      
      ############ KM
      
      n.high <- sum(risk.group, na.rm=T)
      n.low <- sum(!risk.group, na.rm=T)
      
      sdf <- survdiff(Surv(bcr.time, bcr.status) ~ risk.group)
      p.val <- pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE)
      #p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
      
      hr = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
      upper95 = exp(log(hr) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
      lower95 = exp(log(hr) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
      
      stats <- c(stats, c(hr, lower95, upper95, p.val))
      
      statsList[[i]] <- stats
      
    }
  }
}


statsDF <- do.call(rbind, statsList)
statsDF
statsDF <- data.frame(statsDF, stringsAsFactors = F)
statsDF

colnames(statsDF) <- c('Dataset', 'Model', 'Signature', 'C', 'TD.AUC', 'Cox.HR', 'Cox.Lower95',
                       'Cox.Upper95', 'Cox.P', 'KM.HR', 'KM.Lower95', 'KM.Upper95', 'KM.P')

statsDF[,4:13] <- apply(statsDF[,4:13], 2, as.numeric)

dataForBoxPlot <- statsDF

dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxNet'] <- 'CoxLasso'
dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxNetAlpha0'] <- 'CoxRidge'


### SuperPC
dataForBoxPlot$Model[dataForBoxPlot$Model=='SuperPC'] <- 'SuperPC0.3'
dataForBoxPlot$Model <- factor(dataForBoxPlot$Model,
                               levels=c('CoxPH','CoxLasso','CoxRidge','SuperPC0','SuperPC0.1','SuperPC0.3'))

dataForBoxPlot$Dataset <- factor(dataForBoxPlot$Dataset, levels=datasets)
dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, levels=signatures)

saveRDS(dataForBoxPlot, file='report/Summary/Summary_Survival_SuperPC_Intra_Dataset.RDS')


### plsRcox
dataForBoxPlot$Model <- factor(dataForBoxPlot$Model,
                               levels=c('CoxPH','CoxLasso','CoxRidge','plsRcox1','plsRcox2','plsRcox3'))

dataForBoxPlot$Dataset <- factor(dataForBoxPlot$Dataset, levels=datasets)
dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, levels=signatures)

saveRDS(dataForBoxPlot, file='report/Summary/Summary_Survival_plsRcox_Intra_Dataset.RDS')


### RandomForest
dataForBoxPlot$Model <- factor(dataForBoxPlot$Model,
                               levels=c('CoxPH','CoxLasso','CoxRidge','RandomForest10',
                                        'RandomForest20','RandomForest50','RandomForest100'))

dataForBoxPlot$Dataset <- factor(dataForBoxPlot$Dataset, levels=datasets)
dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, levels=signatures)

saveRDS(dataForBoxPlot, file='report/Summary/Summary_Survival_RandomForest_Intra_Dataset.RDS')


med <- dataForBoxPlot %>% group_by(Dataset, Model) %>% 
  summarise(med=median(C, na.rm=T))

med <- dataForBoxPlot %>% group_by(Dataset, Signature) %>% 
  summarise(med=median(C, na.rm=T))

med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(C, na.rm=T))

View(med)

##### Models

ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  geom_boxplot(aes(fill=Model),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  facet_wrap(~Dataset, nrow=2) +
  #ylim(0,5) +
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('C Index')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
