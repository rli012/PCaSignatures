
####################################################################
###              Summary of Survial Analysis Models              ###
####################################################################

library(survivalROC)
library(survival)
library(survminer)
library(survcomp)
library(ggplot2)

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

###### CV10

datasets <- c('TCGA_PRAD','GSE107299','GSE21034','DKFZ2018','GSE54460','GSE70768','GSE70769',
              'GSE94767','E_MTAB_6128','GSE116918_BCR') # , 'GSE116918_Metastasis'


models <- c('CoxPH', 'CoxNetAlpha0', 'SuperPC') #, 'SurvivalSVM'
models <- c('CoxPH', 'CoxNetAlpha0', 'CoxNet', 'SuperPC0', 'SuperPC0.1', 'SuperPC') #, 'SurvivalSVM'

# Published signatures
signatures <- c('Agell','Bibikova','Bismar','Decipher','Ding','Glinsky','Irshad',
                'Jennifer','Jia','Kamoun','Long','Luca','Mo','Nakagawa','Olmos',
                'Oncotype','Penney','Planche','Prolaris','Ramaswamy','Ramos_Montoya',
                'Ross_Adams','Ross_Robert','Sharma','Talantov','Varambally','Wu','Yang',
                'Yu')

i <- 0
statsList <- list()

for (dataset in datasets) {
  
  message(dataset)
  
  for (model in models) {
    
    for (signature.name in signatures) {
      
      message(signature.name)
      
      i <- i + 1
      stats <- c()
      
      # fl.name <- paste0('report/Survival/CV10Scale/Signature/Survival_',
      #                   model, '_', signature.name, '_', dataset, '.csv')
      
      fl.name <- paste0('report/SurvivalC/CV10Scale/Signature/Survival_',
                       model, '_', signature.name, '_', dataset, '.csv')
      
      if (! file.exists(fl.name)) {
        print (paste('NA:', fl.name))
        stats <- c(stats, c(dataset, model, signature.name, rep(NA,10)))
        statsList[[i]] <- stats
        next
      }
      
      res <- read.csv(file=fl.name, stringsAsFactors = F, row.names = 1)
      
      if (sum(res$risk.score==0)>=1/2*length(res$risk.score)) {
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
                       'Cox.Upper96', 'Cox.P', 'KM.HR', 'KM.Lower95', 'KM.Upper95', 'KM.P')

statsDF[,4:13] <- apply(statsDF[,4:13], 2, as.numeric)


dataForBoxPlot <- statsDF

dataForBoxPlot$Dataset <- factor(dataForBoxPlot$Dataset, levels=datasets)
dataForBoxPlot$Model <- factor(dataForBoxPlot$Model, levels=models)
dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, levels=signatures)

# saveRDS(dataForBoxPlot, file='report/Summary/Summary_Survival_CV10.RDS')
# saveRDS(dataForBoxPlot, file='report/Summary/Summary_SurvivalC_CV10.RDS')

library(dplyr)
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
  ylim(0,5) +
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #ylim(0,6) +
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Expression Level (Log'[2]*'CPM)')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
  #theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


###### Signatures

med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(C, na.rm=T))

dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, 
                                   levels=signatures[order(med$med, decreasing = T)])

ggplot(data=dataForBoxPlot, aes(x=Signature, y=C)) +
  geom_boxplot(aes(fill=Signature),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  #facet_wrap(~Dataset, nrow=1) +
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Expression Level (Log'[2]*'CPM)')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
  #theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))




#################################################################

###### Training - Test


datasets <- c('TCGA_PRAD','GSE107299','GSE21034','DKFZ2018','GSE54460','GSE70768','GSE70769',
              'GSE94767','E_MTAB_6128','GSE116918_BCR') # , 'GSE116918_Metastasis'


models <- c('CoxPH', 'CoxNetAlpha0', 'SuperPC') #, 'SurvivalSVM'
models <- c('CoxPH', 'CoxNetAlpha0', 'CoxNet', 'SuperPC0', 'SuperPC0.1', 'SuperPC') #, 'SurvivalSVM'

# Published signatures
signatures <- c('Agell','Bibikova','Bismar','Decipher','Ding','Glinsky','Irshad',
                'Jennifer','Jia','Kamoun','Long','Luca','Mo','Nakagawa','Olmos',
                'Oncotype','Penney','Planche','Prolaris','Ramaswamy','Ramos_Montoya',
                'Ross_Adams','Ross_Robert','Sharma','Talantov','Varambally','Wu','Yang',
                'Yu')

i <- 0
statsList <- list()

for (training.set in datasets) {
  
  message(training.set)
  
  test.sets <- datasets[which(!datasets %in% training.set)]

  for (test.set in test.sets) {
    
    for (model in models) {
      
      for (signature.name in signatures) {
        
        message(signature.name)
        
        i <- i + 1
        stats <- c()
        
        fl.name <- paste0('report/Survival/TrainTestScale/Signature/Survival_',
                          model, '_', signature.name, '_', training.set, '_', test.set, '_TrainTest_Signature.csv')
 
        # fl.name <- paste0('report/SurvivalC/TrainTestScale/Signature/Survival_', 
        #                   model, '_', signature.name, '_', training.set, '_', test.set, '_TrainTest_Signature.csv')
        
        if (! file.exists(fl.name)) {
          print (paste('NA:', fl.name))
          stats <- c(stats, c(training.set, test.set, model, signature.name, rep(NA,10)))
          statsList[[i]] <- stats
          next
        }
        
        res <- read.csv(file=fl.name, stringsAsFactors = F, row.names = 1)
        
        if (sum(res$risk.score==0)>=1/2*length(res$risk.score)) {
          print (paste('Zero:', fl.name))
          stats <- c(stats, c(training.set, test.set, model, signature.name, rep(NA,10)))
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
        
        stats <- c(stats, c(training.set, test.set, model, signature.name, c.index))
        
        
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
}


statsDF <- do.call(rbind, statsList)
statsDF
statsDF <- data.frame(statsDF, stringsAsFactors = F)
statsDF

colnames(statsDF) <- c('Training', 'Test', 'Model', 'Signature', 'C', 'TD.AUC', 'Cox.HR', 'Cox.Lower95',
                       'Cox.Upper96', 'Cox.P', 'KM.HR', 'KM.Lower95', 'KM.Upper95', 'KM.P')

statsDF[,5:14] <- apply(statsDF[,5:14], 2, as.numeric)


dataForBoxPlot <- statsDF

dataForBoxPlot$Training <- factor(dataForBoxPlot$Training, levels=datasets)
dataForBoxPlot$Test <- factor(dataForBoxPlot$Test, levels=datasets)
dataForBoxPlot$Model <- factor(dataForBoxPlot$Model, levels=models)
dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, levels=signatures)

# saveRDS(dataForBoxPlot, file='report/Summary/Summary_Survival_Training_Test.RDS')
# saveRDS(dataForBoxPlot, file='report/Summary/Summary_SurvivalC_Training_Test.RDS')


library(dplyr)
med <- dataForBoxPlot %>% group_by(Training, Model) %>% 
  summarise(med=median(C, na.rm=T))

med <- dataForBoxPlot %>% group_by(Training, Signature) %>% 
  summarise(med=median(C, na.rm=T))

med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(C, na.rm=T))

View(med)


##### Models

ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  geom_boxplot(aes(fill=Model),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  facet_wrap(~Training, nrow=2) +
  #ylim(0,10)+
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Expression Level (Log'[2]*'CPM)')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


###### Signatures

med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(C, na.rm=T))

dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, 
                                   levels=signatures[order(med$med, decreasing = T)])

ggplot(data=dataForBoxPlot, aes(x=Signature, y=C)) +
  geom_boxplot(aes(fill=Signature),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  facet_wrap(~Training, nrow=1) +
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Expression Level (Log'[2]*'CPM)')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))



#################################################################

###### Explore Individual Analysis


datasets <- c('TCGA_PRAD','GSE107299','GSE21034','DKFZ2018','GSE54460','GSE70768','GSE70769',
              'GSE94767','E_MTAB_6128','GSE116918_BCR') # , 'GSE116918_Metastasis'


models <- c('CoxPH', 'CoxNet', 'SuperPC') #, 'SurvivalSVM'
models <- c('CoxPH', 'CoxNetAlpha0', 'CoxNet', 'SuperPC0', 'SuperPC0.1', 'SuperPC') #, 'SurvivalSVM'

# Published signatures
signatures <- c('Agell','Bibikova','Bismar','Decipher','Ding','Glinsky','Irshad',
                'Jennifer','Jia','Kamoun','Long','Luca','Mo','Nakagawa','Olmos',
                'Oncotype','Penney','Planche','Prolaris','Ramaswamy','Ramos_Montoya',
                'Ross_Adams','Ross_Robert','Sharma','Talantov','Varambally','Wu','Yang',
                'Yu')


dataset <- datasets[1]
model <- models[2]
signature.name <- signatures[4]

res <- read.csv(file=paste0('report/Survival/CV10Scale/Signature/Survival_', model, '_', signature.name, '_', dataset, '.csv'),
                stringsAsFactors = F, row.names = 1)

risk.score <- res$risk.score
bcr.time <- res$bcr.time
bcr.status <- res$bcr.status

############ C index

c <- concordance.index(x=risk.score, 
                       surv.time=bcr.time, 
                       surv.event=bcr.status, 
                       #cl=risk.group,
                       method="noether")

c$c.index


############ Time-dependent ROC

cutoff <- 5*12
nobs <- length(risk.score)
roc <- survivalROC(Stime=bcr.time,
                   status=bcr.status,
                   marker = risk.score,
                   method = 'NNE',
                   predict.time = cutoff,
                   span = 0.25*nobs^(-0.20))

roc

FPR <- roc$FP
TPR <- roc$TP

dataForROCPlot <- data.frame(FPR,TPR)

ggplot(dataForROCPlot,aes(x=FPR,y=TPR))+geom_line(size = 1, alpha = 1,color='red')+
  labs(x = "False positive rate (1-specificity)",y = "True positive rate (sensitivity)")+
  #geom_segment(x=0,y=0,xend=1,yend=1, color='black') + 
  geom_abline(slope = 1, intercept = 0) +
  xlim(0,1) + ylim(0,1) +
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black'),
        panel.background = element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) +
  theme(strip.text.x = element_text(size = 12, colour = "black", angle=0)) +
  geom_text(data =NULL, aes(x = 0.85, y = 0.2, label = paste0('AUC = ',round(roc$AUC,3)), group=NULL),
            size=4.5)


############ CoxPH

coxtest <- coxph(Surv(bcr.time, bcr.status) ~ risk.score)
summcph <- summary(coxtest)

coeffs <- c(summcph$coefficients[,2], summcph$conf.int[,3:4], 
            summcph$coefficients[,5])
coeffs


############ KM Plot

risk.group <- res$risk.group
bcr.time <- res$bcr.time
bcr.status <- res$bcr.status

n.high <- sum(risk.group, na.rm=T)
n.low <- sum(!risk.group, na.rm=T)

sdf <- survdiff(Surv(bcr.time, bcr.status) ~ risk.group)
p.val <- pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE)
#p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)

hr = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
upper95 = exp(log(hr) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
lower95 = exp(log(hr) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))

hr <- format(hr, digits = 2, nsmall=2)
upper95 <- format(upper95, digits = 2, nsmall=2)
lower95 <- format(lower95, digits = 2, nsmall=2)

p.val <- ifelse(p.val >= 0.01, formatC(p.val, digits = 2), 
                formatC(p.val, format = "e", digits = 2))

hr
lower95
upper95
p.val


label.hr <- paste('HR = ', hr, ' (', lower95, ' - ', upper95, ')', sep='')
label.p <- paste('P Value = ', p.val, sep='')

km.surv.data <- data.frame(bcr.time, bcr.status, risk.group, stringsAsFactors = F)

fit <- survfit(Surv(bcr.time, bcr.status) ~ risk.group, data=km.surv.data)

lgd.xpos <- 0.27
lgd.ypos = 0.3

p.xpos = max(km.surv.data$bcr.time, na.rm=TRUE)/25
p.ypos = 0.07

#title <- 'PFR10YR'
type <- 'Relapse-free Survival'

plt <- ggsurvplot(fit, data=km.surv.data, pval = paste0(label.hr, '\n', label.p), pval.coord = c(p.xpos, p.ypos),
                  pval.size=5.5,
                  font.main = c(16, 'bold', 'black'), conf.int = FALSE, 
                  #title = title,
                  legend = c(lgd.xpos, lgd.ypos), 
                  #color = c('blue', 'green'),
                  palette= c(google.blue, google.red),
                  legend.labs = c(paste('Low Risk (N=',n.low,')',sep=''), 
                                  paste('High Risk (N=',n.high,')',sep='')),  
                  legend.title='Group',
                  xlab = paste(type,'(months)'), ylab = 'Survival probability',
                  font.x = c(20), font.y = c(20), ylim=c(0,1), #16
                  ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                              panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(),
                                              #panel.border = element_rect(colour='black'),
                                              panel.border = element_blank(),
                                              panel.background = element_blank(),
                                              legend.text = element_text(size=16),#14
                                              legend.title = element_text(size=16),
                                              #axis.title = element_text(size=30),
                                              axis.text = element_text(size=18, color='black')))

print (plt[[1]])




##############################################################

###### Check Files

# Published signatures
signatures <- c('Agell','Bibikova','Bismar','Decipher','Ding','Glinsky','Irshad',
                'Jennifer','Jia','Kamoun','Long','Luca','Mo','Nakagawa','Olmos',
                'Oncotype','Penney','Planche','Prolaris','Ramaswamy','Ramos_Montoya',
                'Ross_Adams','Ross_Robert','Sharma','Talantov','Varambally','Wu','Yang',
                'Yu')

models <- c('CoxPH', 'CoxNet', 'SuperPC') #, 'SurvivalSVM'
models <- c('CoxPH', 'CoxNetAlpha0', 'CoxNet', 'SuperPC0', 'SuperPC0.1', 'SuperPC') #, 'SurvivalSVM'

datasets <- c('TCGA_PRAD','GSE107299','GSE21034','DKFZ2018','GSE54460','GSE70768','GSE70769',
              'GSE94767','E_MTAB_6128','GSE116918_BCR') # , 'GSE116918_Metastasis'

for (dataset in datasets) {
  message(dataset)
  
  for (model in models) {
    for (signature.name in signatures) {
      
      fl <- paste0('report/Survival/CV10Scale/Signature/Survival_', model, '_', signature.name, '_', dataset, '.csv')
      
      if (! file.exists(fl)) {
        print (fl)
        
      }
    }
  }
}


for (training.set in datasets) {
  message(training.set)
  
  test.sets <- datasets[which(!datasets %in% training.set)]
  
  for (test.set in test.sets) {
    
    for (model in models) {
      for (signature.name in signatures) {
        
        fl <- paste0('report/Survival/TrainTestScale/Signature/Survival_', 
                     model, '_', signature.name, '_', training.set, '_', test.set, '_TrainTest_Signature.csv')
        
        if (! file.exists(fl)) {
          print (fl)
          
        }
      }
    }
  }
}

