
####################################################################
###                   Comparison of All Models                   ###
####################################################################

setwd('~/bigdata/PCa/')

library(dplyr)
library(ggplot2)

library(colorspace)
library(RColorBrewer)


#################################################################

###### Intra-Dataset


###### Compare Survival vs. Classification (Same samples)

survc <- readRDS(file='report/Summary/Summary_SurvivalC_Intra_Dataset.RDS')
clsf <- readRDS(file='report/Summary/Summary_Classification_Intra_Dataset.RDS')

survc$Model <- paste0(survc$Model, ' (C)')

dataForBoxPlot <- data.frame(rbind(survc, clsf[,-c(4:9)]), stringsAsFactors = F)

dataForBoxPlot$Model <- factor(dataForBoxPlot$Model, 
                               levels = c('CoxPH (C)','CoxLasso (C)','CoxRidge (C)',
                                          'SuperPC (C)','plsRcox (C)','RandomForest (C)',
                                          'glmnet','svmLinear','svmRadial','svmPoly','rf','pls',
                                          'lda','xgbLinear','xgbTree'))


##### Models

ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  geom_boxplot(aes(fill=Model),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  facet_wrap(~Dataset, nrow=2) +
  geom_vline(xintercept = 6.5, linetype='dashed') +
  #ylim(0,5) +
  #geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #ylim(0,6) +
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Concordance Index')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))




###### Compare Surv vs. SurvC

surv <- readRDS(file='report/Summary/Summary_Survival_Intra_Dataset.RDS')
survc <- readRDS(file='report/Summary/Summary_SurvivalC_Intra_Dataset.RDS')

survc$Model <- paste0(survc$Model, ' (C)')

dataForBoxPlot <- data.frame(rbind(surv, survc), stringsAsFactors = F)

dataForBoxPlot$Model <- factor(dataForBoxPlot$Model, 
                               levels = c('CoxPH','CoxLasso','CoxRidge',
                                          'SuperPC','plsRcox','RandomForest',
                                          'CoxPH (C)','CoxLasso (C)','CoxRidge (C)',
                                          'SuperPC (C)','plsRcox (C)','RandomForest (C)'))


med <- dataForBoxPlot %>% group_by(Dataset, Model) %>% 
  summarise(med=median(C, na.rm=T))
View(med)

##### Models

ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  geom_boxplot(aes(fill=Model),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  facet_wrap(~Dataset, nrow=2) +
  geom_vline(xintercept = 6.5, linetype='dashed') +
  #ylim(0,5) +
  #geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #ylim(0,6) +
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Concordance Index')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))




###### Compare Survival vs. Classification

surv <- readRDS(file='report/Summary/Summary_Survival_Intra_Dataset.RDS')
clsf <- readRDS(file='report/Summary/Summary_Classification_Intra_Dataset.RDS')

dataForBoxPlot <- data.frame(rbind(surv, clsf[,-c(4:9)]), stringsAsFactors = F)

dataForBoxPlot$Model <- factor(dataForBoxPlot$Model, 
                               levels = c('CoxPH','CoxLasso','CoxRidge',
                                          'SuperPC','plsRcox','RandomForest',
                                          'glmnet','svmLinear','svmRadial','svmPoly','rf','pls',
                                          'lda','xgbLinear','xgbTree'))


##### Models

ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  geom_boxplot(aes(fill=Model),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  facet_wrap(~Dataset, nrow=2) +
  geom_vline(xintercept = 6.5, linetype='dashed') +
  #ylim(0,5) +
  #geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #ylim(0,6) +
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Concordance Index')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))



#=======================================================================================#
###### Compare SurvivalC vs. Compare Survival vs. Classification

surv <- readRDS(file='report/Summary/Summary_Survival_Intra_Dataset.RDS')
survc <- readRDS(file='report/Summary/Summary_SurvivalC_Intra_Dataset.RDS')
clsf <- readRDS(file='report/Summary/Summary_Classification_Intra_Dataset.RDS')

survc$Model <- paste0(survc$Model, ' (C)')

dataForBoxPlot <- data.frame(rbind(surv, survc, clsf[,-c(4:9)]), stringsAsFactors = F)

dataForBoxPlot$Model <- factor(dataForBoxPlot$Model, 
                               levels = c('CoxPH','CoxLasso','CoxRidge',
                                          'SuperPC','plsRcox','RandomForest',
                                          'CoxPH (C)','CoxLasso (C)','CoxRidge (C)',
                                          'SuperPC (C)','plsRcox (C)','RandomForest (C)',
                                          'glmnet','svmLinear','svmRadial','svmPoly','rf','pls',
                                          'lda','xgbLinear','xgbTree'))


##### Models

ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  geom_boxplot(aes(fill=Model),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  facet_wrap(~Dataset, nrow=2) +
  geom_vline(xintercept = c(6.5,12.5), linetype='dashed') +
  #ylim(0,5) +
  #geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #ylim(0,6) +
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Concordance Index')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
#=======================================================================================#



###### Compare Survival Models

surv <- readRDS(file='report/Summary/Summary_Survival_Intra_Dataset.RDS')

dataForBoxPlot <- surv

dataForBoxPlot$Model <- factor(dataForBoxPlot$Model, 
                               levels = c('CoxPH','CoxLasso','CoxRidge',
                                          'SuperPC','plsRcox','RandomForest'))


##### Models

ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  geom_boxplot(aes(fill=Model),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  facet_wrap(~Dataset, nrow=2) +
  #geom_vline(xintercept = 6.5, linetype='dashed') +
  #ylim(0,5) +
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #ylim(0,6) +
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Concordance Index')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))



ggplot(data=dataForBoxPlot, aes(x=Model, y=TD.AUC)) +
  geom_boxplot(aes(fill=Model),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  facet_wrap(~Dataset, nrow=2) +
  #geom_vline(xintercept = 6.5, linetype='dashed') +
  #ylim(0,5) +
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #ylim(0,6) +
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Time-dependent AUC')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))



ggplot(data=dataForBoxPlot, aes(x=Model, y=Cox.HR)) +
  geom_boxplot(aes(fill=Model),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  facet_wrap(~Dataset, nrow=2) +
  #geom_vline(xintercept = 6.5, linetype='dashed') +
  #ylim(0,5) +
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #ylim(0,6) +
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Hazard Ratio (Cox)')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))





########## Compare Signatures

surv <- readRDS(file='report/Summary/Summary_Survival_Intra_Dataset.RDS')

dataForBoxPlot <- surv

models <- c('CoxPH','CoxLasso','CoxRidge','SuperPC','plsRcox','RandomForest')

model <- models[3]
dataForBoxPlot <- dataForBoxPlot[dataForBoxPlot$Model==model,]


###### Signatures

### c Index

med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(C, na.rm=T))

dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, 
                                   levels=med$Signature[order(med$med, decreasing = T)])

ggplot(data=dataForBoxPlot, aes(x=Signature, y=C)) +
  geom_boxplot(aes(fill=Signature),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_line(aes(group=Dataset)) +
  #facet_wrap(~Training, nrow=1) +
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Concordance Index')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


### AUC
med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(TD.AUC, na.rm=T))

dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, 
                                   levels=med$Signature[order(med$med, decreasing = T)])

ggplot(data=dataForBoxPlot, aes(x=Signature, y=TD.AUC)) +
  geom_boxplot(aes(fill=Signature),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_line(aes(group=Dataset)) +
  #facet_wrap(~Training, nrow=1) +
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Time-dependent AUC')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


### Cox HR
med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(Cox.HR, na.rm=T))

sort(dataForBoxPlot$Cox.HR, decreasing = T)
#dataForBoxPlot$Cox.HR[dataForBoxPlot$Cox.HR>10] <- 10

dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, 
                                   levels=med$Signature[order(med$med, decreasing = T)])

ggplot(data=dataForBoxPlot, aes(x=Signature, y=Cox.HR)) +
  geom_boxplot(aes(fill=Signature),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_line(aes(group=Dataset)) +
  #facet_wrap(~Training, nrow=1) +
  #ylim(0,10) +
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Hazard Ratio')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


### KM HR
med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(KM.HR, na.rm=T))

dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, 
                                   levels=med$Signature[order(med$med, decreasing = T)])

ggplot(data=dataForBoxPlot, aes(x=Signature, y=KM.HR)) +
  geom_boxplot(aes(fill=Signature),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_line(aes(group=Dataset)) +
  #facet_wrap(~Training, nrow=1) +
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Hazard Ratio')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


x <- cbind(dataForBoxPlot[dataForBoxPlot$Signature=='Li',c(1:6,9)],
           dataForBoxPlot[dataForBoxPlot$Signature=='Penney',c(1:6,9)])

View(x)



##### Forest

surv <- readRDS(file='report/Summary/Summary_Survival_Intra_Dataset.RDS')

dataForForestPlot <- surv

models <- c('CoxPH','CoxLasso','CoxRidge','SuperPC','plsRcox','RandomForest')

model <- models[3]
dataForForestPlot <- dataForForestPlot[dataForForestPlot$Model==model,]

#dataForForestPlot$P <- ifelse(dataForForestPlot$P >= 0.01, formatC(dataForForestPlot$P, digits = 2), 
#                              formatC(dataForForestPlot$P, format = "e", digits = 2))

#dataForForestPlot$P <- paste0('p = ', dataForForestPlot$P)

datasets <- c('TCGA_PRAD','GSE107299','GSE21034','DKFZ2018','GSE54460','GSE70768','GSE70769',
              'GSE94767','E_MTAB_6128','GSE116918_BCR') # , 'GSE116918_Metastasis'

dataset <- datasets[10]
dataForForestPlot <- dataForForestPlot[dataForForestPlot$Dataset==dataset,]


summary(dataForForestPlot$Cox.HR)


### PLOT

o <- order(dataForForestPlot$Cox.HR, decreasing = F)

dataForForestPlot$Signature <- factor(dataForForestPlot$Signature, 
                                   levels=dataForForestPlot$Signature[o])

ggplot(dataForForestPlot, aes(x=Signature, y=Cox.HR)) +
  #geom_segment(aes(y=dataset, x=lower95.coxph, xend=upper95.coxph, yend=dataset), color='black', size=1) +
  #geom_segment(aes(y=6:1-0.1, x=lower95.coxph, xend=lower95.coxph, yend=6:!+0.1), color='black', size=1) +
  geom_errorbar(aes(ymin=Cox.Lower95, ymax=Cox.Upper95),width=0.1, size=0.8, color='black')+ 
  geom_point(color=google.red, size=3, shape=15) + #facet_grid(.~type) +
  #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(-2.9,-3.12,-1.16,1.2,2.58,2.5), label=P, group=NULL),
  #          size=4.4) +
  #ylim(-11,4) +
  geom_hline(yintercept = 0, linetype='dashed') +
  #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(0.35,0.5,0.2,0.45,0.95,0.55), label=p.coxph, group=NULL),
  #          size=4.4) +
  #scale_y_continuous(trans = 'log10',
  #                   breaks = c(0, 1, 2.5,50,250,7500),
  #                   labels = c(0, 1, 2.5,50,250,7500)) +
  coord_flip()+
  #ylim(0,10) +
  xlab('')+ylab(expression('Log (Hazard Ratio)')) +
  facet_wrap(~Dataset, nrow=2) +
  theme_bw()+
  #theme_set(theme_minimal()) #
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = 'right') +
  theme(axis.title=element_text(size=16),
        axis.text = element_text(color='black', size=12),
        axis.text.x = element_text(angle = 0, hjust=0.5),
        strip.text = element_text(size=14)) +
  theme(axis.line = element_line(colour = "black"),
        axis.line.y = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



### KM

cols = brewer.pal(4, "Reds")
cols = colorRampPalette(cols)(10)
col_fun = colorRampPalette(rev(c(cols[10],'white')), space = "Lab")(100)

o <- order(dataForForestPlot$KM.HR, decreasing = F)

dataForForestPlot$Signature <- factor(dataForForestPlot$Signature, 
                                      levels=dataForForestPlot$Signature[o])

ggplot(dataForForestPlot, aes(x=Signature, y=log(KM.HR))) +
  #geom_segment(aes(y=dataset, x=lower95.coxph, xend=upper95.coxph, yend=dataset), color='black', size=1) +
  #geom_segment(aes(y=6:1-0.1, x=lower95.coxph, xend=lower95.coxph, yend=6:!+0.1), color='black', size=1) +
  geom_errorbar(aes(ymin=log(KM.Lower95), ymax=log(KM.Upper95)),width=0.1, size=0.8, color='black')+ 
  #geom_point(color=google.red, size=3, shape=15) + #facet_grid(.~type) +
  geom_point(aes(color=KM.P), size=3, shape=15) + #facet_grid(.~type) +
  #scale_colour_gradientn(limits=c(0,0.05),
  #                       colors= c(google.red, 'white', google.blue)) + #, na.value='grey'
  scale_colour_gradientn(limits=c(0,0.05),
                         colors= c(cols[10], cols[1])) + #, na.value='grey'
  
  #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(-2.9,-3.12,-1.16,1.2,2.58,2.5), label=P, group=NULL),
  #          size=4.4) +
  geom_hline(yintercept = 0, linetype='dashed') +
  #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(0.35,0.5,0.2,0.45,0.95,0.55), label=p.coxph, group=NULL),
  #          size=4.4) +
  #scale_y_continuous(trans = 'log10',
  #                   breaks = c(0, 1, 2.5,50,250,7500),
  #                   labels = c(0, 1, 2.5,50,250,7500)) +
  coord_flip()+
  #ylim(0,10) +
  xlab('')+ylab(expression('Log (Hazard Ratio)')) +
  facet_wrap(~Dataset, nrow=2) +
  theme_bw()+
  #theme_set(theme_minimal()) #
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = 'right') +
  theme(axis.title=element_text(size=16),
        axis.text = element_text(color='black', size=12),
        axis.text.x = element_text(angle = 0, hjust=0.5),
        strip.text = element_text(size=14)) +
  theme(axis.line = element_line(colour = "black"),
        axis.line.y = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())




#################################################################

###### Inter-Dataset


###### Compare Survival vs. Classification

surv <- readRDS(file='report/Summary/Summary_Survival_Inter_Dataset.RDS')
clsf <- readRDS(file='report/Summary/Summary_Classification_Inter_Dataset.RDS')

dataForBoxPlot <- data.frame(rbind(surv, clsf[,-c(5:10)]), stringsAsFactors = F)

dataForBoxPlot$Model <- factor(dataForBoxPlot$Model, 
                               levels = c('CoxPH','CoxLasso','CoxRidge',
                                          'SuperPC','plsRcox','RandomForest',
                                          'glmnet','svmLinear','svmRadial','svmPoly','rf','pls',
                                          'lda','xgbLinear','xgbTree'))

##### Models

ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  geom_boxplot(aes(fill=Model),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  facet_wrap(~Training, nrow=2) +
  geom_vline(xintercept = 6.5, linetype='dashed') +
  #ylim(0,5) +
  #geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #ylim(0,6) +
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Concordance Index')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))



###### Compare Survival Models

surv <- readRDS(file='report/Summary/Summary_Survival_Inter_Dataset.RDS')

dataForBoxPlot <- surv

dataForBoxPlot$Model <- factor(dataForBoxPlot$Model, 
                               levels = c('CoxPH','CoxLasso','CoxRidge',
                                          'SuperPC','plsRcox','RandomForest'))

##### Models

ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  geom_boxplot(aes(fill=Model),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  facet_wrap(~Training, nrow=2) +
  #geom_vline(xintercept = 6.5, linetype='dashed') +
  #ylim(0,5) +
  #geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #ylim(0,6) +
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Concordance Index')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))



ggplot(data=dataForBoxPlot, aes(x=Model, y=TD.AUC)) +
  geom_boxplot(aes(fill=Model),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  facet_wrap(~Training, nrow=2) +
  #geom_vline(xintercept = 6.5, linetype='dashed') +
  #ylim(0,5) +
  #geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #ylim(0,6) +
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Time-dependent AUC')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))



ggplot(data=dataForBoxPlot, aes(x=Model, y=Cox.HR)) +
  geom_boxplot(aes(fill=Model),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  facet_wrap(~Training, nrow=2) +
  #geom_vline(xintercept = 6.5, linetype='dashed') +
  #ylim(0,5) +
  #geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #ylim(0,6) +
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Hazard Ratio (Cox)')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))





########## Compare Signatures

surv <- readRDS(file='report/Summary/Summary_Survival_Inter_Dataset.RDS')

dataForBoxPlot <- surv

models <- c('CoxPH','CoxLasso','CoxRidge','SuperPC','plsRcox','RandomForest')
model <- models[3]
dataForBoxPlot <- dataForBoxPlot[dataForBoxPlot$Model==model,]


###### Signatures

### c Index

med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(C, na.rm=T))

dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, 
                                   levels=med$Signature[order(med$med, decreasing = T)])

ggplot(data=dataForBoxPlot, aes(x=Signature, y=C)) +
  geom_boxplot(aes(fill=Signature),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_line(aes(group=Dataset)) +
  #facet_wrap(~Training, nrow=1) +
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Concordance Index')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


### AUC
med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(TD.AUC, na.rm=T))

dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, 
                                   levels=med$Signature[order(med$med, decreasing = T)])

ggplot(data=dataForBoxPlot, aes(x=Signature, y=TD.AUC)) +
  geom_boxplot(aes(fill=Signature),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_line(aes(group=Dataset)) +
  #facet_wrap(~Training, nrow=1) +
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Time-dependent AUC')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


### Cox HR
med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(Cox.HR, na.rm=T))

dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, 
                                   levels=med$Signature[order(med$med, decreasing = T)])

ggplot(data=dataForBoxPlot, aes(x=Signature, y=Cox.HR)) +
  geom_boxplot(aes(fill=Signature),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_line(aes(group=Dataset)) +
  #facet_wrap(~Training, nrow=1) +
  #ylim(0,10) +
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Hazard Ratio')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


### KM HR
med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(KM.HR, na.rm=T))

dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, 
                                   levels=med$Signature[order(med$med, decreasing = T)])

ggplot(data=dataForBoxPlot, aes(x=Signature, y=KM.HR)) +
  geom_boxplot(aes(fill=Signature),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_line(aes(group=Dataset)) +
  #facet_wrap(~Training, nrow=1) +
  ylim(0,10)+
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Hazard Ratio')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


x <- cbind(dataForBoxPlot[dataForBoxPlot$Signature=='Li',c(1:7,10)],
           dataForBoxPlot[dataForBoxPlot$Signature=='Penney',c(1:7,10)])

View(x)



###### Single dataset

### c Index

surv <- readRDS(file='report/Summary/Summary_Survival_Inter_Dataset.RDS')

dataForBoxPlot <- surv

datasets <- c('TCGA_PRAD','GSE107299','GSE21034','DKFZ2018','GSE54460','GSE70768','GSE70769',
              'GSE94767','E_MTAB_6128','GSE116918_BCR') # , 'GSE116918_Metastasis'


models <- c('CoxPH','CoxLasso','CoxRidge','SuperPC','plsRcox','RandomForest')

model <- models[3]
dataForForestPlot <- dataForForestPlot[dataForForestPlot$Model==model,]

dataset <- datasets[1]
dataForForestPlot <- dataForForestPlot[dataForForestPlot$Training==dataset,]


med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(C, na.rm=T))

dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, 
                                   levels=med$Signature[order(med$med, decreasing = T)])

ggplot(data=dataForBoxPlot, aes(x=Signature, y=C)) +
  geom_boxplot(aes(fill=Signature),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_line(aes(group=Dataset)) +
  facet_wrap(~Training, nrow=2) +
  geom_jitter(size=2, width=0.05, color='black') + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y=expression('Concordance Index')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))




###### Forest Plot

surv <- readRDS(file='report/Summary/Summary_Survival_Inter_Dataset.RDS')

dataForForestPlot <- surv

models <- c('CoxPH','CoxLasso','CoxRidge','SuperPC','plsRcox','RandomForest')

model <- models[3]
dataForForestPlot <- dataForForestPlot[dataForForestPlot$Model==model,]

#dataForForestPlot$P <- ifelse(dataForForestPlot$P >= 0.01, formatC(dataForForestPlot$P, digits = 2), 
#                              formatC(dataForForestPlot$P, format = "e", digits = 2))

#dataForForestPlot$P <- paste0('p = ', dataForForestPlot$P)

datasets <- c('TCGA_PRAD','GSE107299','GSE21034','DKFZ2018','GSE54460','GSE70768','GSE70769',
              'GSE94767','E_MTAB_6128','GSE116918_BCR') # , 'GSE116918_Metastasis'

training.set <- datasets[3]
dataForForestPlot <- dataForForestPlot[dataForForestPlot$Training==training.set,]


test.set <- datasets[5]
dataForForestPlot <- dataForForestPlot[dataForForestPlot$Test==test.set,]


summary(dataForForestPlot$Cox.HR)


### PLOT

o <- order(dataForForestPlot$Test, dataForForestPlot$Cox.HR, decreasing = F)

dataForForestPlot$Signature <- factor(dataForForestPlot$Signature, 
                                      levels=dataForForestPlot$Signature[o])

ggplot(dataForForestPlot, aes(x=Signature, y=Cox.HR)) +
  #geom_segment(aes(y=dataset, x=lower95.coxph, xend=upper95.coxph, yend=dataset), color='black', size=1) +
  #geom_segment(aes(y=6:1-0.1, x=lower95.coxph, xend=lower95.coxph, yend=6:!+0.1), color='black', size=1) +
  geom_errorbar(aes(ymin=Cox.Lower95, ymax=Cox.Upper95),width=0.1, size=0.8, color='black')+ 
  geom_point(color=google.red, size=3, shape=15) + #facet_grid(.~type) +
  #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(-2.9,-3.12,-1.16,1.2,2.58,2.5), label=P, group=NULL),
  #          size=4.4) +
  #ylim(-11,4) +
  geom_hline(yintercept = 0, linetype='dashed') +
  #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(0.35,0.5,0.2,0.45,0.95,0.55), label=p.coxph, group=NULL),
  #          size=4.4) +
  #scale_y_continuous(trans = 'log10',
  #                   breaks = c(0, 1, 2.5,50,250,7500),
  #                   labels = c(0, 1, 2.5,50,250,7500)) +
  coord_flip()+
  #ylim(0,10) +
  xlab('')+ylab(expression('Hazard Ratio')) +
  facet_wrap(~Test, nrow=2, scales = 'free_y') +
  theme_bw()+
  #theme_set(theme_minimal()) #
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = 'right') +
  theme(axis.title=element_text(size=16),
        axis.text = element_text(color='black', size=12),
        axis.text.x = element_text(angle = 0, hjust=0.5),
        strip.text = element_text(size=14)) +
  theme(axis.line = element_line(colour = "black"),
        axis.line.y = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



### KM

cols = brewer.pal(4, "Reds")
cols = colorRampPalette(cols)(10)
col_fun = colorRampPalette(rev(c(cols[10],'white')), space = "Lab")(100)

o <- order(dataForForestPlot$KM.HR, decreasing = F)

dataForForestPlot$Signature <- factor(dataForForestPlot$Signature, 
                                      levels=dataForForestPlot$Signature[o])

ggplot(dataForForestPlot, aes(x=Signature, y=log(KM.HR))) +
  #geom_segment(aes(y=dataset, x=lower95.coxph, xend=upper95.coxph, yend=dataset), color='black', size=1) +
  #geom_segment(aes(y=6:1-0.1, x=lower95.coxph, xend=lower95.coxph, yend=6:!+0.1), color='black', size=1) +
  geom_errorbar(aes(ymin=log(KM.Lower95), ymax=log(KM.Upper95)),width=0.1, size=0.8, color='black')+ 
  #geom_point(color=google.red, size=3, shape=15) + #facet_grid(.~type) +
  geom_point(aes(color=KM.P), size=3, shape=15) + #facet_grid(.~type) +
  #scale_colour_gradientn(limits=c(0,0.05),
  #                       colors= c(google.red, 'white', google.blue)) + #, na.value='grey'
  scale_colour_gradientn(limits=c(0,0.05),
                         colors= c(cols[10], cols[1])) + #, na.value='grey'
  
  #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(-2.9,-3.12,-1.16,1.2,2.58,2.5), label=P, group=NULL),
  #          size=4.4) +
  geom_hline(yintercept = 0, linetype='dashed') +
  #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(0.35,0.5,0.2,0.45,0.95,0.55), label=p.coxph, group=NULL),
  #          size=4.4) +
  #scale_y_continuous(trans = 'log10',
  #                   breaks = c(0, 1, 2.5,50,250,7500),
  #                   labels = c(0, 1, 2.5,50,250,7500)) +
  coord_flip()+
  #ylim(0,10) +
  xlab('')+ylab(expression('Log (Hazard Ratio)')) +
  facet_wrap(~Test, nrow=2) +
  theme_bw()+
  #theme_set(theme_minimal()) #
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = 'right') +
  theme(axis.title=element_text(size=16),
        axis.text = element_text(color='black', size=12),
        axis.text.x = element_text(angle = 0, hjust=0.5),
        strip.text = element_text(size=14)) +
  theme(axis.line = element_line(colour = "black"),
        axis.line.y = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



