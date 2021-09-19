
####################################################################
###                   Comparison of All Models                   ###
####################################################################

setwd('~/Documents/Publications/PCa/')

library(dplyr)
library(ggplot2)

library(colorspace)
library(RColorBrewer)

google.red <- '#ea4235'
google.yellow <- '#fabd03'
google.green <- '#34a853'
google.blue <- '#4286f5'

###########################################################################

###### Compare Survival vs. Classification
## Same sample, more sample
## Intra inter
surv <- readRDS(file='figures/Summary_Survival_Intra_Dataset.RDS')
survc <- readRDS(file='figures/Summary_SurvivalC_Intra_Dataset.RDS')
clsf <- readRDS(file='figures/Summary_Classification_Intra_Dataset.RDS')


## Inter
surv <- readRDS(file='figures/Summary_Survival_Inter_Dataset.RDS')
survc <- readRDS(file='figures/Summary_SurvivalC_Inter_Dataset.RDS')
clsf <- readRDS(file='figures/Summary_Classification_Inter_Dataset.RDS')

survc <- surv

#surv$Type <- 'Survival'
survc$Type <- 'Survival'
clsf$Type <- 'Classification'

View(survc)
View(clsf)

#survc$Model <- paste0(survc$Model, ' (C)')

colnames(survc) == colnames(clsf[,-c(4:9)])

# Intra
dataForBoxPlot <- data.frame(rbind(survc, clsf[,-c(4:9)]), stringsAsFactors = F)
dataForBoxPlot$Dataset <- as.character(dataForBoxPlot$Dataset)

# Inter
dataForBoxPlot <- data.frame(rbind(survc, clsf[,-c(5:10)]), stringsAsFactors = F)
dataForBoxPlot$Dataset <- as.character(dataForBoxPlot$Training)


dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='TCGA_PRAD'] <- 'TCGA-PRAD'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE107299'] <- 'CPC-Gene'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE21034'] <- 'Taylor'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='DKFZ2018'] <- 'DKFZ'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE70768'] <- 'Cambridge'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE70769'] <- 'Stockholm'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE94767'] <- 'CancerMap'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='E_MTAB_6128'] <- 'CIT'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE116918_BCR'] <- 'Belfast'


dataForBoxPlot$Dataset <- factor(dataForBoxPlot$Dataset,
                                 levels = c('TCGA-PRAD','CPC-Gene','Taylor','DKFZ','GSE54460',
                                            'Cambridge','Stockholm','CancerMap','CIT','Belfast'))


unique(dataForBoxPlot$Model)

dataForBoxPlot$Model <- as.character(dataForBoxPlot$Model)

dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxPH'] <- 'CoxPH'
dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxRidge'] <- 'Cox-Ridge'
dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxLasso'] <- 'Cox-Lasso'
dataForBoxPlot$Model[dataForBoxPlot$Model=='SuperPC'] <- 'SuperPC'
dataForBoxPlot$Model[dataForBoxPlot$Model=='plsRcox'] <- 'Cox-PLS'
dataForBoxPlot$Model[dataForBoxPlot$Model=='RandomForest'] <- 'RSF'

dataForBoxPlot$Model[dataForBoxPlot$Model=='glmnet'] <- 'Elastic net'
dataForBoxPlot$Model[dataForBoxPlot$Model=='svmLinear'] <- 'SVM-Linear'
dataForBoxPlot$Model[dataForBoxPlot$Model=='svmPoly'] <- 'SVM-Poly'
dataForBoxPlot$Model[dataForBoxPlot$Model=='svmRadial'] <- 'SVM-RBF'
dataForBoxPlot$Model[dataForBoxPlot$Model=='rf'] <- 'RF'
dataForBoxPlot$Model[dataForBoxPlot$Model=='pls'] <- 'PLS'
dataForBoxPlot$Model[dataForBoxPlot$Model=='lda'] <- 'LDA'
dataForBoxPlot$Model[dataForBoxPlot$Model=='xgbLinear'] <- 'XGBoost-Linear'
dataForBoxPlot$Model[dataForBoxPlot$Model=='xgbTree'] <- 'XGBoost-Tree'

dataForBoxPlot$Model <- factor(dataForBoxPlot$Model, 
                               levels = c('CoxPH','Cox-Ridge','Cox-Lasso',
                                          'SuperPC','Cox-PLS','RSF',
                                          'Elastic net','SVM-Linear','SVM-Poly',
                                          'SVM-RBF','RF','PLS',
                                          'LDA','XGBoost-Linear','XGBoost-Tree'))
dataForBoxPlot$Model

col <- brewer.pal(9, "YlOrBr")[5]

min(dataForBoxPlot$C, na.rm = T)
max(dataForBoxPlot$C, na.rm = T)

View(dataForBoxPlot)

##### Models, dataset facet

ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  geom_boxplot(aes(fill=Type),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  #stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Type)) +
  scale_fill_manual(values = c('Classification'=google.red, 'Survival'=google.blue)) +
  facet_wrap(~Dataset, nrow=2) +
  geom_vline(xintercept = 6.5, linetype='dashed') +
  ylim(0,1) +
  geom_jitter(size=1, width=0.05, color='black', alpha=0.3) + #darkblue
  labs(x='', y='Concordance Index') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))



##### Models, no facet
colnames(dataForBoxPlot)

# C index


p <- ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  # geom_boxplot(aes(color=Type),
  #              outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
  #              outlier.fill = NA) +
  geom_boxplot(aes(color=Type)) +
  stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Type)) +
  scale_color_manual(values = c('Classification'=google.red, 'Survival'=google.blue)) +
  #facet_wrap(~Dataset, nrow=2) +
  #geom_vline(xintercept = 6.5, linetype='dashed') +
  geom_hline(yintercept = 0.5, linetype='dashed') +
  #ylim(0.15,0.85) + # C
  scale_y_continuous(breaks = seq(0.2,0.8,0.2)) +
  #geom_jitter(size=1, width=0.2, color='black', alpha=0.3) + #darkblue
  labs(x='', y='Concordance Index') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col)) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
p

topptx(p, filename = 'figures/Figure_Intra_C_Survival_vs_Classification_Small_Sample_Size.pptx', 
       width=6, height = 3.5)

topptx(p, filename = 'figures/Figure_Inter_C_Survival_vs_Classification_Small_Sample_Size.pptx', 
       width=6, height = 3.5)



# TD AUC
p <- ggplot(data=dataForBoxPlot, aes(x=Model, y=TD.AUC)) +
  # geom_boxplot(aes(color=Type),
  #              outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
  #              outlier.fill = NA) +
  geom_boxplot(aes(color=Type)) +
  stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Type)) +
  scale_color_manual(values = c('Classification'=google.red, 'Survival'=google.blue)) +
  #facet_wrap(~Dataset, nrow=2) +
  #geom_vline(xintercept = 6.5, linetype='dashed') +
  geom_hline(yintercept = 0.5, linetype='dashed') +
  #ylim(0,0.85) + # C
  ylim(0,1) + # C
  #geom_jitter(size=1, width=0.2, color='black', alpha=0.3) + #darkblue
  #labs(x='', y='Concordance Index') +
  labs(x='', y='AUC') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col)) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
p

topptx(p, filename = 'figures/Figure_Intra_AUC_Survival_vs_Classification_Small_Sample_Size.pptx', 
       width=6, height = 3.5)

topptx(p, filename = 'figures/Figure_Inter_AUC_Survival_vs_Classification_Small_Sample_Size.pptx', 
       width=6, height = 3.5)


max(log(dataForBoxPlot$KM.HR), na.rm = T)
min(log(dataForBoxPlot$KM.HR), na.rm = T)

# log(HR)
p <- ggplot(data=dataForBoxPlot, aes(x=Model, y=log(KM.HR))) +
  # geom_boxplot(aes(color=Type),
  #              outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
  #              outlier.fill = NA) +
  geom_boxplot(aes(color=Type)) +
  stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Type)) +
  scale_color_manual(values = c('Classification'=google.red, 'Survival'=google.blue)) +
  #facet_wrap(~Dataset, nrow=2) +
  #geom_vline(xintercept = 6.5, linetype='dashed') +
  geom_hline(yintercept = 0, linetype='dashed') +
  scale_y_continuous(breaks = c(-2,3,1)) +
  #ylim(-3.5,3.5) + # C
  #geom_jitter(size=1, width=0.2, color='black', alpha=0.3) + #darkblue
  #labs(x='', y='Concordance Index') +
  labs(x='', y='Log(Hazard Ratio)') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col)) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
p

topptx(p, filename = 'figures/Figure_Intra_HR_Survival_vs_Classification_Small_Sample_Size.pptx', 
       width=6, height = 3.5)

topptx(p, filename = 'figures/Figure_Inter_HR_Survival_vs_Classification_Small_Sample_Size.pptx', 
       width=6, height = 3.5)

summ <- dataForBoxPlot %>% group_by(Type, Model) %>%
  summarise(Med=median(TD.AUC, na.rm = T))

View(summ)



dataForBoxPlot <- data.frame(Score=c(dataForBoxPlot$C, dataForBoxPlot$TD.AUC, log(dataForBoxPlot$KM.HR)),
                             Model=rep(dataForBoxPlot$Model,3),
                             Type=rep(dataForBoxPlot$Type,3),
                             Metrics=rep(c('Concordance Index','Time-dependent ROC','KM Survival Analysis'), 
                                         each = nrow(dataForBoxPlot)),
                             Dataset=rep(dataForBoxPlot$Dataset, 3))


dataForBoxPlot$Metrics <- factor(dataForBoxPlot$Metrics, 
                                 levels=c('Concordance Index','Time-dependent ROC','KM Survival Analysis'))



# C index, TD AUC, log(HR)
ggplot(data=dataForBoxPlot, aes(x=Model, y=Score)) +
  # geom_boxplot(aes(color=Type),
  #              outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
  #              outlier.fill = NA) +
  geom_boxplot(aes(color=Type)) +
  stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Type)) +
  scale_color_manual(values = c('Classification'=google.red, 'Survival'=google.blue)) +
  facet_wrap(~ Metrics, ncol=1, scales = 'free_y', strip.position = 'right') +
  #geom_vline(xintercept = 6.5, linetype='dashed') +
  #geom_hline(yintercept = c('C'=0.5,'TD.AUC'=0.5,'KM.HR'=0), linetype='dashed') +
  #ylim(-2.2,2.2) + # C
  #geom_jitter(size=1, width=0.2, color='black', alpha=0.3) + #darkblue
  #labs(x='', y='Concordance Index') +
  labs(x='', y='log(Hazard Ratio)') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col)) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))



###########################################################################

###### Compare Survival vs. Survival C
## Intra inter

surv <- readRDS(file='figures/Summary_Survival_Intra_Dataset.RDS')
survc <- readRDS(file='figures/Summary_SurvivalC_Intra_Dataset.RDS')
#clsf <- readRDS(file='figures/Summary_Classification_Intra_Dataset.RDS')

surv <- readRDS(file='figures/Summary_Survival_Inter_Dataset.RDS')
survc <- readRDS(file='figures/Summary_SurvivalC_Inter_Dataset.RDS')
#clsf <- readRDS(file='figures/Summary_Classification_Inter_Dataset.RDS')

survc$Type <- 'Partial'
surv$Type <- 'Complete'

View(survc)
View(clsf)

#survc$Model <- paste0(survc$Model, ' (C)')


# Intra
dataForBoxPlot <- data.frame(rbind(surv, survc), stringsAsFactors = F)
dataForBoxPlot$Dataset <- as.character(dataForBoxPlot$Dataset)

# Inter
dataForBoxPlot <- data.frame(rbind(surv, survc), stringsAsFactors = F)
dataForBoxPlot$Dataset <- as.character(dataForBoxPlot$Training)


dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='TCGA_PRAD'] <- 'TCGA-PRAD'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE107299'] <- 'CPC-Gene'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE21034'] <- 'Taylor'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='DKFZ2018'] <- 'DKFZ'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE70768'] <- 'Cambridge'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE70769'] <- 'Stockholm'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE94767'] <- 'CancerMap'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='E_MTAB_6128'] <- 'CIT'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE116918_BCR'] <- 'Belfast'


dataForBoxPlot$Dataset <- factor(dataForBoxPlot$Dataset,
                                 levels = c('TCGA-PRAD','CPC-Gene','Taylor','DKFZ','GSE54460',
                                            'Cambridge','Stockholm','CancerMap','CIT','Belfast'))

dataForBoxPlot$Model <- as.character(dataForBoxPlot$Model)

dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxPH'] <- 'CoxPH'
dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxRidge'] <- 'Cox-Ridge'
dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxLasso'] <- 'Cox-Lasso'
dataForBoxPlot$Model[dataForBoxPlot$Model=='SuperPC'] <- 'SuperPC'
dataForBoxPlot$Model[dataForBoxPlot$Model=='plsRcox'] <- 'Cox-PLS'
dataForBoxPlot$Model[dataForBoxPlot$Model=='RandomForest'] <- 'RSF'

dataForBoxPlot$Model <- factor(dataForBoxPlot$Model, 
                               levels = c('CoxPH','Cox-Ridge','Cox-Lasso',
                                          'SuperPC','Cox-PLS','RSF'))
dataForBoxPlot$Model


col <- brewer.pal(9, "YlOrBr")[5]


min(dataForBoxPlot$C, na.rm = T)
max(dataForBoxPlot$C, na.rm = T)

##### Models, dataset facet

ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  geom_boxplot(aes(fill=Type),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  #stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Type)) +
  scale_fill_manual(values = c('Classification'=google.red, 'Survival'=google.blue)) +
  facet_wrap(~Dataset, nrow=2) +
  geom_vline(xintercept = 6.5, linetype='dashed') +
  ylim(0,1) +
  geom_jitter(size=1, width=0.05, color='black', alpha=0.3) + #darkblue
  labs(x='', y='Concordance Index') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))




##### Models, no facet
colnames(dataForBoxPlot)
View(dataForBoxPlot)

pos <- position_dodge(0.85)

# C index

p <- ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  # geom_boxplot(aes(color=Type),
  #              outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
  #              outlier.fill = NA) +
  geom_boxplot(aes(color=Type), position = pos) +
  stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Type),
               position = pos) +
  scale_color_manual(values = c('Partial'=google.yellow, 'Complete'=google.blue)) +
  #facet_wrap(~Dataset, nrow=2) +
  #geom_vline(xintercept = 6.5, linetype='dashed') +
  geom_hline(yintercept = 0.5, linetype='dashed') +
  #ylim(0.15,0.85) + # C
  ylim(0,1) + # C
  #geom_jitter(size=1, width=0.2, color='black', alpha=0.3) + #darkblue
  labs(x='', y='Concordance Index') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col)) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
p

topptx(p, filename = 'figures/Figure_Intra_C_Survival_vs_SurvivalC.pptx', 
       width=6, height = 3.5)

topptx(p, filename = 'figures/Figure_Inter_C_Survival_vs_SurvivalC.pptx', 
       width=6, height = 3.5)


# TD AUC
p <- ggplot(data=dataForBoxPlot, aes(x=Model, y=TD.AUC)) +
  # geom_boxplot(aes(color=Type),
  #              outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
  #              outlier.fill = NA) +
  geom_boxplot(aes(color=Type), position = pos) +
  stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Type),
               position = pos) +
  scale_color_manual(values = c('Partial'=google.yellow, 'Complete'=google.blue)) +
  #facet_wrap(~Dataset, nrow=2) +
  #geom_vline(xintercept = 6.5, linetype='dashed') +
  geom_hline(yintercept = 0.5, linetype='dashed') +
  #ylim(0,0.85) + # C
  ylim(0,1) + # C
  #geom_jitter(size=1, width=0.2, color='black', alpha=0.3) + #darkblue
  #labs(x='', y='Concordance Index') +
  labs(x='', y='AUC') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col)) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
p

topptx(p, filename = 'figures/Figure_Intra_AUC_Survival_vs_SurvivalC.pptx', 
       width=6, height = 3.5)

topptx(p, filename = 'figures/Figure_Inter_AUC_Survival_vs_SurvivalC.pptx', 
       width=6, height = 3.5)


max(log(dataForBoxPlot$KM.HR), na.rm = T)
min(log(dataForBoxPlot$KM.HR), na.rm = T)

# log(HR)
p <- ggplot(data=dataForBoxPlot, aes(x=Model, y=log(KM.HR))) +
  # geom_boxplot(aes(color=Type),
  #              outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
  #              outlier.fill = NA) +
  geom_boxplot(aes(color=Type), position = pos) +
  stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Type),
               position = pos) +
  scale_color_manual(values = c('Partial'=google.yellow, 'Complete'=google.blue)) +
  #facet_wrap(~Dataset, nrow=2) +
  #geom_vline(xintercept = 6.5, linetype='dashed') +
  geom_hline(yintercept = 0, linetype='dashed') +
  #scale_y_continuous(breaks = seq(-3,3,1)) +
  #ylim(-3.5,3.5) + # C
  #geom_jitter(size=1, width=0.2, color='black', alpha=0.3) + #darkblue
  #labs(x='', y='Concordance Index') +
  labs(x='', y='log(Hazard Ratio)') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col)) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
p

topptx(p, filename = 'figures/Figure_Intra_HR_Survival_vs_SurvivalC.pptx', 
       width=6, height = 3.5)

topptx(p, filename = 'figures/Figure_Inter_HR_Survival_vs_SurvivalC.pptx', 
       width=6, height = 3.5)

summ <- dataForBoxPlot %>% group_by(Type, Model) %>%
  summarise(Med=median(C, na.rm = T))

View(summ)



dataForBoxPlot <- data.frame(Score=c(dataForBoxPlot$C, dataForBoxPlot$TD.AUC, log(dataForBoxPlot$KM.HR)),
                             Model=rep(dataForBoxPlot$Model,3),
                             Type=rep(dataForBoxPlot$Type,3),
                             Metrics=rep(c('Concordance Index','Time-dependent ROC','KM Survival Analysis'), 
                                         each = nrow(dataForBoxPlot)),
                             Dataset=rep(dataForBoxPlot$Dataset, 3))


dataForBoxPlot$Metrics <- factor(dataForBoxPlot$Metrics, 
                                 levels=c('Concordance Index','Time-dependent ROC','KM Survival Analysis'))



# C index, TD AUC, log(HR)
ggplot(data=dataForBoxPlot, aes(x=Model, y=Score)) +
  # geom_boxplot(aes(color=Type),
  #              outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
  #              outlier.fill = NA) +
  geom_boxplot(aes(color=Type)) +
  stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Type)) +
  scale_color_manual(values = c('Classification'=google.red, 'Survival'=google.blue)) +
  facet_wrap(~ Metrics, ncol=1, scales = 'free_y', strip.position = 'right') +
  #geom_vline(xintercept = 6.5, linetype='dashed') +
  #geom_hline(yintercept = c('C'=0.5,'TD.AUC'=0.5,'KM.HR'=0), linetype='dashed') +
  #ylim(-2.2,2.2) + # C
  #geom_jitter(size=1, width=0.2, color='black', alpha=0.3) + #darkblue
  #labs(x='', y='Concordance Index') +
  labs(x='', y='log(Hazard Ratio)') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col)) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))





###########################################################################

###### Compare Survival Models in each dataset
## Intra inter


# Intra
surv <- readRDS(file='figures/Summary_Survival_Intra_Dataset.RDS')
dataForBoxPlot <- surv
dataForBoxPlot$Dataset <- as.character(dataForBoxPlot$Dataset)

# Inter
surv <- readRDS(file='figures/Summary_Survival_Inter_Dataset.RDS')
dataForBoxPlot <- surv
dataForBoxPlot$Dataset <- as.character(dataForBoxPlot$Training)


dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='TCGA_PRAD'] <- 'TCGA-PRAD'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE107299'] <- 'CPC-Gene'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE21034'] <- 'Taylor'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='DKFZ2018'] <- 'DKFZ'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE70768'] <- 'Cambridge'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE70769'] <- 'Stockholm'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE94767'] <- 'CancerMap'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='E_MTAB_6128'] <- 'CIT'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE116918_BCR'] <- 'Belfast'


dataForBoxPlot$Dataset <- factor(dataForBoxPlot$Dataset,
                                 levels = c('TCGA-PRAD','CPC-Gene','Taylor','DKFZ','GSE54460',
                                            'Cambridge','Stockholm','CancerMap','CIT','Belfast'))

dataForBoxPlot$Model <- as.character(dataForBoxPlot$Model)

dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxPH'] <- 'CoxPH'
dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxRidge'] <- 'Cox-Ridge'
dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxLasso'] <- 'Cox-Lasso'
dataForBoxPlot$Model[dataForBoxPlot$Model=='SuperPC'] <- 'SuperPC'
dataForBoxPlot$Model[dataForBoxPlot$Model=='plsRcox'] <- 'Cox-PLS'
dataForBoxPlot$Model[dataForBoxPlot$Model=='RandomForest'] <- 'RSF'

dataForBoxPlot$Model <- factor(dataForBoxPlot$Model, 
                               levels = c('CoxPH','Cox-Ridge','Cox-Lasso',
                                          'SuperPC','Cox-PLS','RSF'))
dataForBoxPlot$Model

col <- brewer.pal(9, "YlOrBr")[4]


min(dataForBoxPlot$C, na.rm = T)
max(dataForBoxPlot$C, na.rm = T)

##### Models, dataset facet, C

p <- ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  geom_boxplot(aes(color=Model),
               #outlier.shape = NA, 
               outlier.size = 1.2, #outlier.colour = 'black',
               #outlier.fill = NA
  ) +
  stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Model)) +
  #scale_fill_manual(values = c('Classification'=google.red, 'Survival'=google.blue)) +
  facet_wrap(~Dataset, nrow=2) +
  geom_hline(yintercept = 0.5, linetype='dashed') +
  ylim(0,1) +
  #geom_jitter(size=0.2, width=0.15, color='black', alpha=0.3) + #darkblue
  labs(x='', y='Concordance Index') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=12, face='bold')) +
  theme(strip.background = element_rect(fill=col, size = 0.5, color=col)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))

p
topptx(p, filename = 'figures/Figure_Intra_C_Survival.pptx', 
       width=8.5, height = 5)

topptx(p, filename = 'figures/Figure_Inter_C_Survival.pptx', 
       width=8.5, height = 5)




##### Models, dataset facet, AUC

p <- ggplot(data=dataForBoxPlot, aes(x=Model, y=TD.AUC)) +
  geom_boxplot(aes(color=Model),
               #outlier.shape = NA, 
               outlier.size = 1.2, #outlier.colour = 'black',
               #outlier.fill = NA
  ) +
  stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Model)) +
  #scale_fill_manual(values = c('Classification'=google.red, 'Survival'=google.blue)) +
  facet_wrap(~Dataset, nrow=2) +
  geom_hline(yintercept = 0.5, linetype='dashed') +
  ylim(0,1) +
  #geom_jitter(size=1, width=0.1, color='black', alpha=0.3) + #darkblue
  labs(x='', y='AUC') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=12, face='bold')) +
  theme(strip.background = element_rect(fill=col, size = 0.5, color=col)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))

p

topptx(p, filename = 'figures/Figure_Intra_AUC_Survival.pptx', 
       width=8.5, height = 5)

topptx(p, filename = 'figures/Figure_Inter_AUC_Survival.pptx', 
       width=8.5, height = 5)

##### Models, dataset facet, log(HR)

p <- ggplot(data=dataForBoxPlot, aes(x=Model, y=log(KM.HR))) +
  geom_boxplot(aes(color=Model),
               #outlier.shape = NA, 
               outlier.size = 1.2, #outlier.colour = 'black',
               #outlier.fill = NA
  ) +
  stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Model)) +
  #scale_fill_manual(values = c('Classification'=google.red, 'Survival'=google.blue)) +
  facet_wrap(~Dataset, nrow=2) +
  geom_hline(yintercept = 0, linetype='dashed') +
  #ylim(0,1) +
  #geom_jitter(size=1, width=0.1, color='black', alpha=0.3) + #darkblue
  labs(x='', y='Log(Hazard Ratio)') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=12, face='bold')) +
  theme(strip.background = element_rect(fill=col, size = 0.5, color=col)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))

p

topptx(p, filename = 'figures/Figure_Intra_HR_Survival.pptx', 
       width=8.5, height = 5)

topptx(p, filename = 'figures/Figure_Inter_HR_Survival.pptx', 
       width=8.5, height = 5)

library(dplyr)

summ <- dataForBoxPlot %>% group_by(Dataset, Model) %>%
  summarise(Med=median(TD.AUC, na.rm = T))

View(summ)
View(dataForBoxPlot)

write.table(summ, file='figures/Summary_Survival_Inter_Dataset_C_Index.txt',
            sep = '\t', quote = F, row.names = F)


##### Models, dataset facet, other

ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  geom_violin(aes(fill=Model),
              outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
              outlier.fill = NA) +
  geom_boxplot(fill='white', width = 0.3,
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  #stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Type)) +
  #scale_fill_manual(values = c('Classification'=google.red, 'Survival'=google.blue)) +
  facet_wrap(~Dataset, nrow=2) +
  geom_hline(yintercept = 0.5, linetype='dashed') +
  ylim(0,1) +
  #geom_jitter(size=1, width=0.1, color='black', alpha=0.3) + #darkblue
  labs(x='', y='Concordance Index') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col, size = 0.5, color=col)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))





ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  geom_boxplot(aes(fill=Model),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  #stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Type)) +
  #scale_fill_manual(values = c('Classification'=google.red, 'Survival'=google.blue)) +
  facet_wrap(~Dataset, nrow=2) +
  geom_hline(yintercept = 0.5, linetype='dashed') +
  ylim(0,1) +
  #geom_jitter(size=1, width=0.15, color='black', alpha=0.3) + #darkblue
  labs(x='', y='Concordance Index') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col, size = 0.5, color=col)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))



ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  geom_boxplot(aes(color=Model),
               outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA) +
  stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Model)) +
  #scale_fill_manual(values = c('Classification'=google.red, 'Survival'=google.blue)) +
  facet_wrap(~Dataset, nrow=2) +
  geom_hline(yintercept = 0.5, linetype='dashed') +
  ylim(0,1) +
  geom_jitter(size=1, width=0.15, color='black', alpha=0.3) + #darkblue
  labs(x='', y='Concordance Index') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col, size = 0.5, color=col)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))



#########################################################################


########## Compare Signatures

surv <- readRDS(file='figures/Summary_Survival_Intra_Dataset.RDS')
surv <- readRDS(file='figures/Summary_Survival_Inter_Dataset.RDS')

tr <- readRDS(file='figures/Summary_Survival_Transcriptome_Intra_Dataset.RDS')
tr <- readRDS(file='figures/Summary_Survival_Transcriptome_Inter_Dataset.RDS')


dataForBoxPlot <- data.frame(rbind(surv, tr), stringsAsFactors = F)

dataForBoxPlot <- surv

dataForBoxPlot$Model <- as.character(dataForBoxPlot$Model)

dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxPH'] <- 'CoxPH'
dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxRidge'] <- 'Cox-Ridge'
dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxLasso'] <- 'Cox-Lasso'
dataForBoxPlot$Model[dataForBoxPlot$Model=='SuperPC'] <- 'SuperPC'
dataForBoxPlot$Model[dataForBoxPlot$Model=='plsRcox'] <- 'Cox-PLS'
dataForBoxPlot$Model[dataForBoxPlot$Model=='RandomForest'] <- 'RSF'


models <- c('CoxPH','Cox-Ridge','Cox-Lasso','SuperPC','Cox-PLS','RSF')

model <- models[2]
dataForBoxPlot <- dataForBoxPlot[dataForBoxPlot$Model==model,]

dataForBoxPlot$Signature <- as.character(dataForBoxPlot$Signature)

dataForBoxPlot$Signature[dataForBoxPlot$Signature=='Prolaris'] <- 'Cuzick'
dataForBoxPlot$Signature[dataForBoxPlot$Signature=='Decipher'] <- 'Erho'
dataForBoxPlot$Signature[dataForBoxPlot$Signature=='Oncotype'] <- 'Klein'
dataForBoxPlot$Signature[dataForBoxPlot$Signature=='Jennifer'] <- 'Sinnott'
dataForBoxPlot$Signature[dataForBoxPlot$Signature=='Ramos_Montoya'] <- 'Ramos-Montoya'
dataForBoxPlot$Signature[dataForBoxPlot$Signature=='Ross_Adams'] <- 'Ross-Adams'
dataForBoxPlot$Signature[dataForBoxPlot$Signature=='Ross_Robert'] <- 'Ross'


sort(unique(as.character(dataForBoxPlot$Signature)))

###### Signatures

### c Index


med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(C, na.rm=T))

dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, 
                                   levels=med$Signature[order(med$med, decreasing = T)])

cols <- ifelse(med$Signature[order(med$med, decreasing = T)]=='Transcriptome',google.red, 'black')


p <- ggplot(data=dataForBoxPlot, aes(x=Signature, y=C)) +
  geom_boxplot(aes(fill=Signature),
               outlier.shape = NA, 
               #outlier.size = 2, 
               #outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_line(aes(group=Dataset)) +
  geom_hline(yintercept = 0.5, linetype='dashed') +
  #facet_wrap(~Training, nrow=1) +
  #ylim(0,1) +
  geom_jitter(size=1.5, width=0.1, color='black', alpha = 0.3) + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y='Concordance Index') +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = cols),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col, size = 0.5, color=col)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))

p

topptx(p, filename = 'figures/Figure_Intra_C_Signature_CoxRidge.pptx', 
       width=8, height = 4)

topptx(p, filename = 'figures/Figure_Inter_C_Signature_CoxRidge.pptx', 
       width=8, height = 4)

topptx(p, filename = 'figures/Figure_Intra_C_Signature_CoxPLS.pptx', 
       width=8, height = 4)

topptx(p, filename = 'figures/Figure_Inter_C_Signature_CoxPLS.pptx', 
       width=8, height = 4)

topptx(p, filename = 'figures/Figure_Intra_C_Signature_CoxRidge_Transcriptome.pptx', 
       width=8, height = 4)


topptx(p, filename = 'figures/Figure_Inter_C_Signature_CoxRidge_Transcriptome.pptx', 
       width=8, height = 4)


### AUC
med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(TD.AUC, na.rm=T))

dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, 
                                   levels=med$Signature[order(med$med, decreasing = T)])

cols <- ifelse(med$Signature[order(med$med, decreasing = T)]=='Transcriptome',google.red, 'black')


p <- ggplot(data=dataForBoxPlot, aes(x=Signature, y=TD.AUC)) +
  geom_boxplot(aes(fill=Signature),
               outlier.shape = NA, 
               #outlier.size = 2, 
               #outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_line(aes(group=Dataset)) +
  geom_hline(yintercept = 0.5, linetype='dashed') +
  #facet_wrap(~Training, nrow=1) +
  #ylim(0,1) +
  geom_jitter(size=1.5, width=0.1, color='black', alpha = 0.3) + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y='AUC') +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = cols),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col, size = 0.5, color=col)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
p

topptx(p, filename = 'figures/Figure_Intra_AUC_Signature_CoxRidge.pptx', 
       width=8, height = 4)

topptx(p, filename = 'figures/Figure_Inter_AUC_Signature_CoxRidge.pptx', 
       width=8, height = 4)

topptx(p, filename = 'figures/Figure_Intra_AUC_Signature_CoxPLS.pptx', 
       width=8, height = 4)

topptx(p, filename = 'figures/Figure_Inter_AUC_Signature_CoxPLS.pptx', 
       width=8, height = 4)

topptx(p, filename = 'figures/Figure_Intra_AUC_Signature_CoxRidge_Transcriptome.pptx', 
       width=8, height = 4)

topptx(p, filename = 'figures/Figure_Inter_AUC_Signature_CoxRidge_Transcriptome.pptx', 
       width=8, height = 4)

### KM HR
med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(KM.HR, na.rm=T))

dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, 
                                   levels=med$Signature[order(med$med, decreasing = T)])

cols <- ifelse(med$Signature[order(med$med, decreasing = T)]=='Transcriptome',google.red, 'black')

p <- ggplot(data=dataForBoxPlot, aes(x=Signature, y=log(KM.HR))) +
  geom_boxplot(aes(fill=Signature),
               outlier.shape = NA, 
               #outlier.size = 2, 
               #outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_line(aes(group=Dataset)) +
  geom_hline(yintercept = 0, linetype='dashed') +
  #facet_wrap(~Training, nrow=1) +
  #ylim(0,1) +
  geom_jitter(size=1.5, width=0.1, color='black', alpha = 0.3) + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y='Log(Hazard Ratio)') +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = cols),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col, size = 0.5, color=col)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))

p

topptx(p, filename = 'figures/Figure_Intra_KM_Signature_CoxRidge.pptx', 
       width=8, height = 4)

topptx(p, filename = 'figures/Figure_Inter_KM_Signature_CoxRidge.pptx', 
       width=8, height = 4)

topptx(p, filename = 'figures/Figure_Intra_KM_Signature_CoxPLS.pptx', 
       width=8, height = 4)

topptx(p, filename = 'figures/Figure_Inter_KM_Signature_CoxPLS.pptx', 
       width=8, height = 4)

topptx(p, filename = 'figures/Figure_Intra_KM_Signature_CoxRidge_Transcriptome.pptx', 
       width=8, height = 4)

topptx(p, filename = 'figures/Figure_Inter_KM_Signature_CoxRidge_Transcriptome.pptx', 
       width=8, height = 4)

### Cox HR (Not included in the manuscript)
med <- dataForBoxPlot %>% group_by(Signature) %>% 
  summarise(med=median(Cox.HR, na.rm=T))

sort(dataForBoxPlot$Cox.HR, decreasing = T)
#dataForBoxPlot$Cox.HR[dataForBoxPlot$Cox.HR>10] <- 10

dataForBoxPlot$Signature <- factor(dataForBoxPlot$Signature, 
                                   levels=med$Signature[order(med$med, decreasing = T)])

ggplot(data=dataForBoxPlot, aes(x=Signature, y=log(Cox.HR))) +
  geom_boxplot(aes(fill=Signature),
               outlier.shape = NA, 
               #outlier.size = 2, 
               #outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_line(aes(group=Dataset)) +
  geom_hline(yintercept = 0, linetype='dashed') +
  #facet_wrap(~Training, nrow=1) +
  ylim(-5,5) +
  geom_jitter(size=1.5, width=0.1, color='black', alpha = 0.3) + #darkblue
  #scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
  labs(x='', y='log(Hazard Ratio)') +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col, size = 0.5, color=col)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


dataForBoxPlot[log(dataForBoxPlot$Cox.HR) < -5,]


x <- cbind(dataForBoxPlot[dataForBoxPlot$Signature=='Li',c(1:6,9)],
           dataForBoxPlot[dataForBoxPlot$Signature=='Penney',c(1:6,9)])

View(x)



dim(dataForForestPlot)
View(dataForForestPlot)

##### Forest

surv <- readRDS(file='figures/Summary_Survival_Intra_Dataset.RDS')

dataForForestPlot <- surv

dataForForestPlot$Model <- as.character(dataForForestPlot$Model)

dataForForestPlot$Model[dataForForestPlot$Model=='CoxPH'] <- 'CoxPH'
dataForForestPlot$Model[dataForForestPlot$Model=='CoxRidge'] <- 'Cox-Ridge'
dataForForestPlot$Model[dataForForestPlot$Model=='CoxLasso'] <- 'Cox-Lasso'
dataForForestPlot$Model[dataForForestPlot$Model=='SuperPC'] <- 'SuperPC'
dataForForestPlot$Model[dataForForestPlot$Model=='plsRcox'] <- 'Cox-PLS'
dataForForestPlot$Model[dataForForestPlot$Model=='RandomForest'] <- 'RSF'

models <- c('CoxPH','Cox-Ridge','Cox-Lasso','SuperPC','Cox-PLS','RSF')


model <- models[2]
dataForForestPlot <- dataForForestPlot[dataForForestPlot$Model==model,]

#dataForForestPlot$P <- ifelse(dataForForestPlot$P >= 0.01, formatC(dataForForestPlot$P, digits = 2), 
#                              formatC(dataForForestPlot$P, format = "e", digits = 2))

#dataForForestPlot$P <- paste0('p = ', dataForForestPlot$P)


dataForForestPlot$Signature <- as.character(dataForForestPlot$Signature)

dataForForestPlot$Signature[dataForForestPlot$Signature=='Prolaris'] <- 'Cuzick'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Decipher'] <- 'Erho'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Oncotype'] <- 'Klein'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Jennifer'] <- 'Sinnott'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Ramos_Montoya'] <- 'Ramos-Montoya'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Ross_Adams'] <- 'Ross-Adams'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Ross_Robert'] <- 'Ross'

# Intra
dataForForestPlot$Dataset <- as.character(dataForForestPlot$Dataset)

# Inter
dataForForestPlot$Dataset <- as.character(dataForForestPlot$Training)

dataForForestPlot$Dataset[dataForForestPlot$Dataset=='TCGA_PRAD'] <- 'TCGA-PRAD'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE107299'] <- 'CPC-Gene'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE21034'] <- 'Taylor'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='DKFZ2018'] <- 'DKFZ'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE70768'] <- 'Cambridge'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE70769'] <- 'Stockholm'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE94767'] <- 'CancerMap'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='E_MTAB_6128'] <- 'CIT'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE116918_BCR'] <- 'Belfast'


dataForForestPlot$Dataset <- factor(dataForForestPlot$Dataset,
                                    levels = c('TCGA-PRAD','CPC-Gene','Taylor','DKFZ','GSE54460',
                                               'Cambridge','Stockholm','CancerMap','CIT','Belfast'))


datasets <- c('TCGA-PRAD','CPC-Gene','Taylor','DKFZ','GSE54460',
              'Cambridge','Stockholm','CancerMap','CIT','Belfast')

dataset <- datasets[1]
dataForForestPlot <- dataForForestPlot[dataForForestPlot$Dataset==dataset,]

dataForForestPlot$Dataset
dataset

summary(dataForForestPlot$Cox.HR)


### PLOT

# o <- order(dataForForestPlot$Cox.HR, decreasing = F)
# 
# dataForForestPlot$Signature <- factor(dataForForestPlot$Signature, 
#                                    levels=dataForForestPlot$Signature[o])
# 
# ggplot(dataForForestPlot, aes(x=Signature, y=Cox.HR)) +
#   #geom_segment(aes(y=dataset, x=lower95.coxph, xend=upper95.coxph, yend=dataset), color='black', size=1) +
#   #geom_segment(aes(y=6:1-0.1, x=lower95.coxph, xend=lower95.coxph, yend=6:!+0.1), color='black', size=1) +
#   geom_errorbar(aes(ymin=Cox.Lower95, ymax=Cox.Upper95),width=0.1, size=0.8, color='black')+ 
#   geom_point(color=google.red, size=3, shape=15) + #facet_grid(.~type) +
#   #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(-2.9,-3.12,-1.16,1.2,2.58,2.5), label=P, group=NULL),
#   #          size=4.4) +
#   #ylim(-11,4) +
#   geom_hline(yintercept = 0, linetype='dashed') +
#   #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(0.35,0.5,0.2,0.45,0.95,0.55), label=p.coxph, group=NULL),
#   #          size=4.4) +
#   #scale_y_continuous(trans = 'log10',
#   #                   breaks = c(0, 1, 2.5,50,250,7500),
#   #                   labels = c(0, 1, 2.5,50,250,7500)) +
#   coord_flip()+
#   #ylim(0,10) +
#   xlab('')+ylab(expression('Log (Hazard Ratio)')) +
#   facet_wrap(~Dataset, nrow=2) +
#   theme_bw()+
#   #theme_set(theme_minimal()) #
#   theme(legend.title = element_blank(),
#         legend.text = element_text(size=14),
#         legend.position = 'right') +
#   theme(axis.title=element_text(size=16),
#         axis.text = element_text(color='black', size=12),
#         axis.text.x = element_text(angle = 0, hjust=0.5),
#         strip.text = element_text(size=14)) +
#   theme(axis.line = element_line(colour = "black"),
#         axis.line.y = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())



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
  geom_errorbar(aes(ymin=log(KM.Lower95), ymax=log(KM.Upper95)),width=0.4, size=0.8, color='black')+ 
  #geom_point(color=google.red, size=3, shape=15) + #facet_grid(.~type) +
  geom_point(color=google.red, size=3, shape=15) + #facet_grid(.~type) +
  #scale_colour_gradientn(limits=c(0,0.05),
  #                       colors= c(google.red, 'white', google.blue)) + #, na.value='grey'
  # scale_colour_gradientn(limits=c(0,0.05),
  #                        colors= c(cols[10], cols[1])) + #, na.value='grey'
  
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
  xlab('')+ylab('Log(Hazard Ratio)') +
  facet_wrap(~Dataset, nrow=2) +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.title=element_text(size=14, face = 'bold'),
        axis.line = element_line(colour = "black"),
        axis.line.y = element_blank(),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col))



####

col <- brewer.pal(9, "YlOrBr")[4]

o <- order(dataForForestPlot$Dataset, dataForForestPlot$KM.HR, decreasing = F)
dataForForestPlot <- dataForForestPlot[o,]

dataForForestPlot$Test <- dataForForestPlot$Dataset

p <- dataForForestPlot %>% 
  mutate(Signature = reorder(Signature, KM.HR)) %>%
  group_by(Test, Signature) %>% 
  arrange(desc(KM.HR)) %>% 
  ungroup() %>% 
  mutate(Signature = factor(paste(Signature, Test, sep = "__"), 
                            levels = rev(paste(Signature, Test, sep = "__")))) %>%
  ggplot(aes(Signature, log(KM.HR))) +
  
  #geom_segment(aes(y=dataset, x=lower95.coxph, xend=upper95.coxph, yend=dataset), color='black', size=1) +
  #geom_segment(aes(y=6:1-0.1, x=lower95.coxph, xend=lower95.coxph, yend=6:!+0.1), color='black', size=1) +
  geom_errorbar(aes(ymin=log(KM.Lower95), ymax=log(KM.Upper95)),width=0.3, size=0.8, color='black')+ 
  #geom_point(color=google.red, size=3, shape=15) + #facet_grid(.~type) +
  geom_point(color=google.red, size=1.5, shape=15) + #facet_grid(.~type) +
  scale_x_discrete(labels = function(x) gsub("__.+$", "", x)) +
  #scale_colour_gradientn(limits=c(0,0.05),
  #                       colors= c(google.red, 'white', google.blue)) + #, na.value='grey'
  # scale_colour_gradientn(limits=c(0,0.05),
  #                        colors= c(cols[10], cols[1])) + #, na.value='grey'
  
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
  xlab('')+ylab('Log(Hazard Ratio)') +
  facet_wrap(~Test, nrow=2, scales = 'free_y') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=10, face = 'bold', color = 'black'),
        axis.title=element_text(size=12, face = 'bold'),
        axis.line = element_line(colour = "black"),
        axis.line.y = element_blank(),
        strip.text = element_text(size=11, face='bold')) +
  theme(strip.background = element_rect(fill=col))

p

topptx(p, filename = 'figures/Figure_Intra_KM_Signature_Each_Dataset.pptx', 
       width=15, height = 9)








##### bubble, Intra

surv <- readRDS(file='figures/Summary_Survival_Intra_Dataset.RDS')

dataForForestPlot <- surv

models <- c('CoxPH','CoxLasso','CoxRidge','SuperPC','plsRcox','RandomForest')

model <- models[5]
dataForForestPlot <- dataForForestPlot[dataForForestPlot$Model==model,]

#dataForForestPlot$P <- ifelse(dataForForestPlot$P >= 0.01, formatC(dataForForestPlot$P, digits = 2), 
#                              formatC(dataForForestPlot$P, format = "e", digits = 2))

#dataForForestPlot$P <- paste0('p = ', dataForForestPlot$P)


dataForForestPlot$Signature <- as.character(dataForForestPlot$Signature)

dataForForestPlot$Signature[dataForForestPlot$Signature=='Prolaris'] <- 'Cuzick'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Decipher'] <- 'Erho'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Oncotype'] <- 'Klein'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Jennifer'] <- 'Sinnott'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Ramos_Montoya'] <- 'Ramos-Montoya'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Ross_Adams'] <- 'Ross-Adams'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Ross_Robert'] <- 'Ross'

# Intra
dataForForestPlot$Dataset <- as.character(dataForForestPlot$Dataset)

dataForForestPlot$Dataset[dataForForestPlot$Dataset=='TCGA_PRAD'] <- 'TCGA-PRAD'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE107299'] <- 'CPC-Gene'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE21034'] <- 'Taylor'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='DKFZ2018'] <- 'DKFZ'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE70768'] <- 'Cambridge'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE70769'] <- 'Stockholm'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE94767'] <- 'CancerMap'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='E_MTAB_6128'] <- 'CIT'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE116918_BCR'] <- 'Belfast'


dataForForestPlot$Dataset <- factor(dataForForestPlot$Dataset,
                                    levels = c('TCGA-PRAD','CPC-Gene','Taylor','DKFZ','GSE54460',
                                               'Cambridge','Stockholm','CancerMap','CIT','Belfast'))

# o <- order(dataForForestPlot$KM.P[dataForForestPlot$Dataset=='TCGA-PRAD'], decreasing = T)
# o
# 
# dataForForestPlot$Signature <- factor(dataForForestPlot$Signature, 
#                                       levels = dataForForestPlot$Signature[dataForForestPlot$Dataset=='TCGA-PRAD'][o])


x <- dataForForestPlot %>% group_by(Signature) %>%
  summarize(Count=sum(KM.P<0.05))

x[order(x$Count, decreasing = T),]

o <- x[order(x$Count, decreasing = T),]$Signature
o

dataForForestPlot$Signature <- factor(dataForForestPlot$Signature, 
                                      levels = rev(o))


dataForBubblePlot <- dataForForestPlot


p <- ggplot(dataForBubblePlot, mapping=aes(x=Signature, y=Dataset, #y=-log10(Benjamini), #y=Fold.Enrichment
                                      color=KM.P,size=log(KM.HR))) +
  geom_point()+ coord_flip() +
  #scale_x_discrete(limits=rev(unique(dataForBubblePlot$Signature))) +
  #scale_x_discrete(limits=Order)+
  scale_colour_gradientn(limits=c(0,0.05),
                         colors= c(google.red, google.blue)) + #
  #facet_wrap(~Comparison) +
  #facet_grid(Regulation~Comparison) + # scales=free
  xlab('')+ylab('') + #ggtitle("") + 
  guides(size = guide_legend(order=2, title = 'Log(HR)'),
         colour = guide_colourbar(order=1, title = 'P Value')) + #'P Value\n(Benjamini)'))
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour='black'),
                   panel.background = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size=20)) +
  theme(axis.text=element_text(size=12, color='black', face = 'bold'),
        axis.text.x =element_text(size=13, color='black', face = 'bold', angle=45, hjust=1),
        axis.title=element_text(size=14, face = 'bold')) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = 'bold')) +
  theme(legend.key.size = unit(0.8,'cm')) +
  theme(strip.background = element_rect(fill=col),
        strip.text = element_text(size=11, face='bold'))

p
topptx(p, filename = 'figures/Figure_Inter_KM_Signature_Each_Dataset_1.pptx', 
       width=5.2, height = 6.5)


#################################################################################


##### bubble, Inter

surv <- readRDS(file='figures/Summary_Survival_Inter_Dataset.RDS')

dataForForestPlot <- surv

dataForForestPlot$Model <- as.character(dataForForestPlot$Model)

dataForForestPlot$Model[dataForForestPlot$Model=='CoxPH'] <- 'CoxPH'
dataForForestPlot$Model[dataForForestPlot$Model=='CoxRidge'] <- 'Cox-Ridge'
dataForForestPlot$Model[dataForForestPlot$Model=='CoxLasso'] <- 'Cox-Lasso'
dataForForestPlot$Model[dataForForestPlot$Model=='SuperPC'] <- 'SuperPC'
dataForForestPlot$Model[dataForForestPlot$Model=='plsRcox'] <- 'Cox-PLS'
dataForForestPlot$Model[dataForForestPlot$Model=='RandomForest'] <- 'RSF'

models <- c('CoxPH','Cox-Ridge','Cox-Lasso','SuperPC','Cox-PLS','RSF')

model <- models[2]
dataForForestPlot <- dataForForestPlot[dataForForestPlot$Model==model,]

max.size <- log(max(dataForForestPlot$KM.HR))
max.size

#dataForForestPlot$P <- ifelse(dataForForestPlot$P >= 0.01, formatC(dataForForestPlot$P, digits = 2), 
#                              formatC(dataForForestPlot$P, format = "e", digits = 2))

#dataForForestPlot$P <- paste0('p = ', dataForForestPlot$P)


dataForForestPlot$Signature <- as.character(dataForForestPlot$Signature)

dataForForestPlot$Signature[dataForForestPlot$Signature=='Prolaris'] <- 'Cuzick'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Decipher'] <- 'Erho'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Oncotype'] <- 'Klein'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Jennifer'] <- 'Sinnott'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Ramos_Montoya'] <- 'Ramos-Montoya'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Ross_Adams'] <- 'Ross-Adams'
dataForForestPlot$Signature[dataForForestPlot$Signature=='Ross_Robert'] <- 'Ross'

training <- c('TCGA_PRAD','GSE107299','GSE21034','DKFZ2018','GSE54460',
              'GSE70768','GSE70769','GSE94767','E_MTAB_6128','GSE116918_BCR')

i <- 10
dataForForestPlot <- dataForForestPlot[dataForForestPlot$Training==training[i],]


# Inter
dataForForestPlot$Dataset <- as.character(dataForForestPlot$Test)

dataForForestPlot$Dataset[dataForForestPlot$Dataset=='TCGA_PRAD'] <- 'TCGA-PRAD'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE107299'] <- 'CPC-Gene'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE21034'] <- 'Taylor'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='DKFZ2018'] <- 'DKFZ'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE70768'] <- 'Cambridge'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE70769'] <- 'Stockholm'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE94767'] <- 'CancerMap'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='E_MTAB_6128'] <- 'CIT'
dataForForestPlot$Dataset[dataForForestPlot$Dataset=='GSE116918_BCR'] <- 'Belfast'


dataForForestPlot$Dataset <- factor(dataForForestPlot$Dataset,
                                    levels = c('TCGA-PRAD','CPC-Gene','Taylor','DKFZ','GSE54460',
                                               'Cambridge','Stockholm','CancerMap','CIT','Belfast'))

# o <- order(dataForForestPlot$KM.P[dataForForestPlot$Dataset=='TCGA-PRAD'], decreasing = T)
# o
# 
# dataForForestPlot$Signature <- factor(dataForForestPlot$Signature, 
#                                       levels = dataForForestPlot$Signature[dataForForestPlot$Dataset=='TCGA-PRAD'][o])


x <- dataForForestPlot %>% group_by(Signature) %>%
  summarize(Count=sum(KM.P<0.05 & KM.HR>1))

#x[order(x$Count, decreasing = T),]

o <- x[order(x$Count, decreasing = T),]$Signature
o

dataForForestPlot$Signature <- factor(dataForForestPlot$Signature, 
                                      levels = rev(o))

dataForBubblePlot <- dataForForestPlot
#View(dataForBubblePlot)

filter <- which(dataForForestPlot$KM.P>0.05 | dataForForestPlot$KM.HR<1)
dataForBubblePlot$KM.HR[filter] <- NA


p <- ggplot(dataForBubblePlot, mapping=aes(x=Signature, y=Dataset, #y=-log10(Benjamini), #y=Fold.Enrichment
                                           color=KM.P,size=log(KM.HR))) +
  geom_point(na.rm = TRUE)+ coord_flip() +
  #scale_x_discrete(limits=rev(unique(dataForBubblePlot$Signature))) +
  #scale_x_discrete(limits=Order)+
  scale_colour_gradientn(limits=c(0,0.05),
                         colors= c(google.red, google.blue)) + #
  scale_size_continuous(limits = c(0,3.1)) +
  #facet_wrap(~Comparison) +
  #facet_grid(Regulation~Comparison) + # scales=free
  xlab('')+ylab('') + #ggtitle("") + 
  guides(size = guide_legend(order=2, title = 'Log(HR)'),
         colour = guide_colourbar(order=1, title = 'P Value')) + #'P Value\n(Benjamini)'))
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour='black'),
                   panel.background = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size=20)) +
  theme(axis.text=element_text(size=12, color='black', face = 'bold'),
        axis.text.x =element_text(size=13, color='black', face = 'bold', angle=45, hjust=1),
        axis.title=element_text(size=14, face = 'bold')) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = 'bold')) +
  theme(legend.key.size = unit(0.6,'cm')) #+
  theme(legend.position = 'none')

p

i <- 'legend'
topptx(p, filename = paste0('figures/Figure_Inter_KM_Signature_Each_Dataset_', i, '.pptx'), 
       width=4, height = 6.5)


filter <- which(dataForForestPlot$KM.P>0.05)
dataForBubblePlot <- dataForForestPlot[-filter,]

dataForBubblePlot <- dataForForestPlot
View(dataForBubblePlot)

o <- order(dataForBubblePlot$Training, dataForBubblePlot$Test, dataForBubblePlot$KM.HR, decreasing = F)
dataForBubblePlot <- dataForBubblePlot[o,]

#dataForBubblePlot$Test <- dataForBubblePlot$Dataset

View(p)

p <- dataForBubblePlot %>% 
  mutate(Signature = reorder(Signature, KM.HR)) %>%
  group_by(Training, Test, Signature) %>% 
  arrange(desc(KM.HR)) %>% 
  ungroup() %>% 
  mutate(Signature = factor(paste(Signature, Training, Test, sep = "__"), 
                            levels = rev(paste(Signature,Training, Test, sep = "__")))) %>%
  
  ggplot(mapping=aes(x=Signature, y=Test, #y=-log10(Benjamini), #y=Fold.Enrichment
                     color=KM.P,size=log(KM.HR))) +
  geom_point()+ coord_flip() + 
  scale_x_discrete(labels = function(x) gsub("__.+$", "", x)) +
  #scale_x_discrete(limits=rev(unique(dataForBubblePlot$Signature))) +
  #scale_x_discrete(limits=Order)+
  scale_colour_gradientn(limits=c(0,0.05), na.value = 'transparent',
                         colors= c(google.red, google.blue)) + #
  facet_wrap(~Training, nrow=2, scales = 'free') +
  #facet_grid(Regulation~Comparison) + # scales=free
  xlab('')+ylab('') + #ggtitle("") + 
  #scale_size_continuous(range = c(0,max.size)) +
  guides(size = guide_legend(order=2, title = 'Log(HR)'),
         colour = guide_colourbar(order=1, title = 'P Value')) + #'P Value\n(Benjamini)'))
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour='black'),
                   panel.background = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size=20)) +
  theme(axis.text=element_text(size=10, color='black', face = 'bold'),
        axis.text.x =element_text(size=12, color='black', face = 'bold', angle=45, hjust=1),
        axis.title=element_text(size=11, face = 'bold')) +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = 'bold')) +
  theme(legend.key.size = unit(0.8,'cm'))

p


topptx(p, filename = 'figures/Figure_Intra_KM_Signature_Each_Dataset.pptx', 
       width=15, height = 9)









#########################################################################################

### Transcriptome

surv <- readRDS(file='figures/Summary_Survival_Intra_Dataset.RDS')
tr <- readRDS(file='figures/Summary_Survival_Transcriptome_Intra_Dataset.RDS')

surv <- readRDS(file='figures/Summary_Survival_Inter_Dataset.RDS')
tr <- readRDS(file='figures/Summary_Survival_Transcriptome_Inter_Dataset.RDS')


dataForBoxPlot <- data.frame(rbind(surv, tr), stringsAsFactors = F)

#models <- c('CoxPH','CoxLasso','CoxRidge','SuperPC','plsRcox','RandomForest')

# model <- models[3]
# dataForBoxPlot <- dataForBoxPlot[dataForBoxPlot$Model==model,]

dataForBoxPlot$Signature <- as.character(dataForBoxPlot$Signature)

dataForBoxPlot$Signature[dataForBoxPlot$Signature=='Prolaris'] <- 'Cuzick'
dataForBoxPlot$Signature[dataForBoxPlot$Signature=='Decipher'] <- 'Erho'
dataForBoxPlot$Signature[dataForBoxPlot$Signature=='Oncotype'] <- 'Klein'
dataForBoxPlot$Signature[dataForBoxPlot$Signature=='Jennifer'] <- 'Sinnott'
dataForBoxPlot$Signature[dataForBoxPlot$Signature=='Ramos_Montoya'] <- 'Ramos-Montoya'
dataForBoxPlot$Signature[dataForBoxPlot$Signature=='Ross_Adams'] <- 'Ross-Adams'
dataForBoxPlot$Signature[dataForBoxPlot$Signature=='Ross_Robert'] <- 'Ross'

dataForBoxPlot$Dataset <- as.character(dataForBoxPlot$Dataset)

dataForBoxPlot$Dataset <- as.character(dataForBoxPlot$Training)

dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='TCGA_PRAD'] <- 'TCGA-PRAD'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE107299'] <- 'CPC-Gene'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE21034'] <- 'Taylor'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='DKFZ2018'] <- 'DKFZ'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE70768'] <- 'Cambridge'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE70769'] <- 'Stockholm'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE94767'] <- 'CancerMap'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='E_MTAB_6128'] <- 'CIT'
dataForBoxPlot$Dataset[dataForBoxPlot$Dataset=='GSE116918_BCR'] <- 'Belfast'

dataForBoxPlot$Test <- as.character(dataForBoxPlot$Test)

dataForBoxPlot$Test[dataForBoxPlot$Test=='TCGA_PRAD'] <- 'TCGA-PRAD'
dataForBoxPlot$Test[dataForBoxPlot$Test=='GSE107299'] <- 'CPC-Gene'
dataForBoxPlot$Test[dataForBoxPlot$Test=='GSE21034'] <- 'Taylor'
dataForBoxPlot$Test[dataForBoxPlot$Test=='DKFZ2018'] <- 'DKFZ'
dataForBoxPlot$Test[dataForBoxPlot$Test=='GSE70768'] <- 'Cambridge'
dataForBoxPlot$Test[dataForBoxPlot$Test=='GSE70769'] <- 'Stockholm'
dataForBoxPlot$Test[dataForBoxPlot$Test=='GSE94767'] <- 'CancerMap'
dataForBoxPlot$Test[dataForBoxPlot$Test=='E_MTAB_6128'] <- 'CIT'
dataForBoxPlot$Test[dataForBoxPlot$Test=='GSE116918_BCR'] <- 'Belfast'


dataForBoxPlot$Dataset <- factor(dataForBoxPlot$Dataset,
                                 levels = c('TCGA-PRAD','CPC-Gene','Taylor','DKFZ','GSE54460',
                                            'Cambridge','Stockholm','CancerMap','CIT','Belfast'))

dataForBoxPlot$Test <- factor(dataForBoxPlot$Test,
                              levels = c('TCGA-PRAD','CPC-Gene','Taylor','DKFZ','GSE54460',
                                         'Cambridge','Stockholm','CancerMap','CIT','Belfast'))


dataForBoxPlot$Model <- as.character(dataForBoxPlot$Model)

dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxPH'] <- 'CoxPH'
dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxRidge'] <- 'Cox-Ridge'
dataForBoxPlot$Model[dataForBoxPlot$Model=='CoxLasso'] <- 'Cox-Lasso'
dataForBoxPlot$Model[dataForBoxPlot$Model=='SuperPC'] <- 'SuperPC'
dataForBoxPlot$Model[dataForBoxPlot$Model=='plsRcox'] <- 'Cox-PLS'
dataForBoxPlot$Model[dataForBoxPlot$Model=='RandomForest'] <- 'RSF'

dataForBoxPlot$Model <- factor(dataForBoxPlot$Model, 
                               levels = c('CoxPH','Cox-Ridge','Cox-Lasso',
                                          'SuperPC','Cox-PLS','RSF'))

sort(unique(as.character(dataForBoxPlot$Signature)))

model <- 'Transcriptome'
dataForBoxPlot <- dataForBoxPlot[dataForBoxPlot$Signature==model,]

View(dataForBoxPlot)


###### Model

### c Index


col <- brewer.pal(9, "YlOrBr")[4]

min(dataForBoxPlot$C, na.rm = T)
max(dataForBoxPlot$C, na.rm = T)

##### Models, dataset facet, C

p <- ggplot(data=dataForBoxPlot, aes(x=Model, y=C)) +
  geom_boxplot(aes(color=Model),
               #outlier.shape = NA, 
               outlier.size = 1.2, #outlier.colour = 'black',
               #outlier.fill = NA
  ) +
  stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Model)) +
  #scale_fill_manual(values = c('Classification'=google.red, 'Survival'=google.blue)) +
  #facet_wrap(~Dataset, nrow=2) +
  geom_hline(yintercept = 0.5, linetype='dashed') +
  ylim(0,1) +
  geom_jitter(size=2, width=0.1, color='black', alpha=0.3) + #darkblue
  labs(x='', y='Concordance Index') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col, size = 0.5, color=col)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))

p

topptx(p, filename = 'figures/Figure_Intra_C_Survival_Model_Comparison_Transcriptome.pptx', 
       width=6, height = 3.5)

topptx(p, filename = 'figures/Figure_Inter_C_Survival_Model_Comparison_Transcriptome.pptx', 
       width=6, height = 3.5)


##### Models, dataset facet, AUC

p <- ggplot(data=dataForBoxPlot, aes(x=Model, y=TD.AUC)) +
  geom_boxplot(aes(color=Model),
               outlier.shape = NA, 
               outlier.size = NA, #outlier.colour = 'black',
               outlier.fill = NA
  ) +
  stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Model)) +
  #scale_fill_manual(values = c('Classification'=google.red, 'Survival'=google.blue)) +
  #facet_wrap(~Dataset, nrow=2) +
  geom_hline(yintercept = 0.5, linetype='dashed') +
  ylim(0,1) +
  geom_jitter(size=2, width=0.1, color='black', alpha=0.3) + #darkblue
  labs(x='', y='AUC') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col, size = 0.5, color=col)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))

p

topptx(p, filename = 'figures/Figure_Intra_AUC_Survival_Model_Comparison_Transcriptome.pptx', 
       width=6, height = 3.5)

topptx(p, filename = 'figures/Figure_Inter_AUC_Survival_Model_Comparison_Transcriptome.pptx', 
       width=6, height = 3.5)


##### Models, dataset facet, log(HR)

p <- ggplot(data=dataForBoxPlot, aes(x=Model, y=log(KM.HR))) +
  geom_boxplot(aes(color=Model),
               outlier.shape = NA, 
               outlier.size = 1.2, #outlier.colour = 'black',
               outlier.fill = NA
  ) +
  stat_summary(fun.y='median', geom="point", shape=23, size=3, aes(color=Model)) +
  #scale_fill_manual(values = c('Classification'=google.red, 'Survival'=google.blue)) +
  #facet_wrap(~Dataset, nrow=2) +
  geom_hline(yintercept = 0, linetype='dashed') +
  #ylim(0,1) +
  geom_jitter(size=2, width=0.1, color='black', alpha=0.3) + #darkblue
  labs(x='', y='Log(Hazard Ratio)') +
  theme_bw()+theme(#axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour='black'),
    panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12, face = 'bold', color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14, face = 'bold'),
        strip.text = element_text(size=14, face='bold')) +
  theme(strip.background = element_rect(fill=col, size = 0.5, color=col)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))

p

topptx(p, filename = 'figures/Figure_Intra_KM_Survival_Model_Comparison_Transcriptome.pptx', 
       width=6, height = 3.5)

topptx(p, filename = 'figures/Figure_Inter_KM_Survival_Model_Comparison_Transcriptome.pptx', 
       width=6, height = 3.5)




summ <- dataForBoxPlot %>% group_by(Dataset, Model) %>%
  summarise(Med=median(C, na.rm = T))

View(summ)

write.table(summ, file='figures/Summary_Survival_Inter_Dataset_C_Index_Transcriptome.txt',
            sep = '\t', quote = F, row.names = F)


######### Stat



datasets <- c('TCGA-PRAD','CPC-Gene','Taylor','DKFZ','GSE54460',
              'Cambridge','Stockholm','CancerMap','CIT','Belfast')

stats <- c()

for (model in models[2:6]) {
  
  print (model)
  
  dataForStats <- dataForBoxPlot[dataForBoxPlot$Model==model,]
  
  for (dataset in datasets) {
    print (dataset)
    
    dataForStats2 <- dataForStats[dataForStats$Dataset==dataset,]
    
    rank.c <- which(dataForStats2$Signature[order(dataForStats2$C, decreasing = T)]=='Transcriptome')
    rank.auc <- which(dataForStats2$Signature[order(dataForStats2$TD.AUC, decreasing = T)]=='Transcriptome')
    rank.hr <- which(dataForStats2$Signature[order(dataForStats2$KM.HR, decreasing = T)]=='Transcriptome')
    
    stats <- rbind(stats, c(model, dataset, rank.c, rank.auc, rank.hr))
    
  }
  
}

colnames(stats) <- c('Model','Dataset','Rank.C','Rank.AUC','Rank.KM.HR')

write.table(stats, file='figures/Comparison_Transcriptome_Rank_Intra_Dataset.txt',
            sep = '\t', quote = F, row.names = F)




idx <- which(apply(dataForBoxPlot, 1, function(x) sum(is.na(x)))>0)
dataForBoxPlot <- dataForBoxPlot[-idx,]


datasets <- c('TCGA-PRAD','CPC-Gene','Taylor','DKFZ','GSE54460',
              'Cambridge','Stockholm','CancerMap','CIT','Belfast')

stats <- c()


for (model in models[2:6]) {
  
  #print (model)
  
  dataForStats <- dataForBoxPlot[dataForBoxPlot$Model==model,]
  
  for (dataset in datasets) {
    #print (dataset)
    
    #dataForStats2 <- dataForStats[dataForStats$Dataset==dataset,]
    
    for (test in datasets) {
      
      dataForStats2 <- dataForStats[which(dataForStats$Dataset==dataset & dataForStats$Test==test),]
      
      rank.c <- which(dataForStats2$Signature[order(dataForStats2$C, decreasing = T)]=='Transcriptome')
      rank.auc <- which(dataForStats2$Signature[order(dataForStats2$TD.AUC, decreasing = T)]=='Transcriptome')
      rank.hr <- which(dataForStats2$Signature[order(dataForStats2$KM.HR, decreasing = T)]=='Transcriptome')
      
      if (length(rank.c)==0) {
        print ('yes')
        next
        
      }
      
      stats <- rbind(stats, c(model, dataset, test, rank.c, rank.auc, rank.hr))
      
    }
    
  }
  
}

View(dataForStats)
View(stats)

colnames(stats) <- c('Model','Training','Test','Rank.C','Rank.AUC','Rank.KM.HR')

write.table(stats, file='figures/Comparison_Transcriptome_Rank_Inter_Dataset.txt',
            sep = '\t', quote = F, row.names = F)


#####################################################

### Rank

library(ggforce)

surv <- readRDS(file='figures/Summary_Survival_Intra_Dataset.RDS')

#dataForParallelPlot <- surv[which(surv$Dataset=='TCGA_PRAD'),]
dataForParallelPlot <- surv

dataForScatterPlot <- dataForParallelPlot %>% group_by(Dataset, Model) %>%
  mutate(Rank=rank(C*-1))

View(dataForScatterPlot)

# data <- reshape2::melt(Titanic)
# data <- gather_set_data(data, 1:4)
# 
# ggplot(data, aes(x, id = id, value = value)) +
#   geom_parallel_sets(aes(fill = Sex), alpha = 0.3, axis.width = 0.1) +
#   geom_parallel_sets_axes(axis.width = 0.1) +
#   geom_parallel_sets_labels(colour = 'white')

dataForScatterPlot <- dataForScatterPlot[dataForScatterPlot$Model %in% c('CoxRidge','plsRcox'),]

dataForScatterPlot

ggplot(data=dataForScatterPlot, aes(x=Model, y=Rank, group=Signature)) +
  geom_line(aes(color=Signature), size=2, alpha=0.3) +
  geom_point(size=2, color=google.red) +
  #scale_color_manual(values=c("#4285F4", "#FBBC05")) +
  #facet_wrap(~groupCol) +
  labs(x='', y=expression('Expression Level (Log'[2]*'Intensity)')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2, color=NULL, group=NULL)) +
  # geom_text(data =anno, aes(x, y, label=label, group=NULL, colour=NULL),
  #           size=4) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 0, hjust = 0),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) +
  #stat_compare_means(comparisons = my_comparisons,
  #                   label = 'p.signif', # p.format
  #                   hide.ns = F,
  #                   #paired = T,
  #                   tip.length = 0.02) + ### NOT IN USE
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))

View(dataForScatterPlot)


#####################################################

### Rank heatmap

datasets <- c('TCGA-PRAD','CPC-Gene','Taylor','DKFZ','GSE54460',
              'Cambridge','Stockholm','CancerMap','CIT','Belfast')

rank.intra <- read.table(file='figures/Comparison_Transcriptome_Rank_Intra_Dataset.txt',
                         sep = '\t', header = T, stringsAsFactors = F)

rank.intra

library(tibble)
rank.intra <- add_column(rank.intra, rank.intra$Dataset, .after = 2)

colnames(rank.intra)[2:3] <- c('Training','Test')


rank.inter <- read.table(file='figures/Comparison_Transcriptome_Rank_Inter_Dataset.txt',
                         sep = '\t', header = T, stringsAsFactors = F)


rank.stats <- rbind(rank.intra, rank.inter)

View(rank.stats)

library(reshape2)

dataForHeatmap <- rank.stats[rank.stats$Model=='CoxRidge',]
View(dataForHeatmap)

dataForHeatmap <- dcast(dataForHeatmap, Training ~ Test, value.var="Rank.C")
dataForHeatmap <- dcast(dataForHeatmap, Training ~ Test, value.var="Rank.AUC")
dataForHeatmap <- dcast(dataForHeatmap, Training ~ Test, value.var="Rank.KM.HR")

dataForHeatmap <- data.frame(dataForHeatmap, stringsAsFactors = F)

rownames(dataForHeatmap) <- dataForHeatmap$Training
dataForHeatmap <- dataForHeatmap[,-1]

dataForHeatmap

colnames(dataForHeatmap) <- gsub('.', '-', colnames(dataForHeatmap), fixed = T)


dataForHeatmap <- dataForHeatmap[datasets,datasets]

annoMatrix <- dataForHeatmap

cols = brewer.pal(4, "Reds")
cols = colorRampPalette(cols)(10)
cols

col_fun = colorRampPalette(rev(c(cols[10],cols[1])), space = "Lab")(30)


ht <- Heatmap(as.matrix(dataForHeatmap),
              #name = 'Expression',
              
              # COLOR
              #col = colorRampPalette(rev(c("red",'white','blue')), space = "Lab")(100),
              col=rev(col_fun),
              na_col = 'gray',
              rect_gp = gpar(col = "grey", lwd = 0.5),
              
              # MAIN PANEL
              column_title = NULL,
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              show_column_dend = FALSE,
              show_row_names = TRUE,
              row_names_side = "left",
              column_names_side = 'top',
              show_column_names = TRUE,
              column_names_rot = 45,
              column_names_gp = gpar(fontsize = 12, fontface='bold'),
              row_names_gp = gpar(fontsize = 12, fontface='bold'),
              #column_names_max_height = unit(3, 'cm'),
              #column_split = factor(phenoData$Day,
              #                      levels=str_sort(unique(phenoData$Day), numeric = T)),
              
              #column_order = rownames(phenoData),
              
              # ANNOTATION
              #top_annotation = topAnnotation,
              
              # ADD TEXT
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(annoMatrix[i, j], x, y, gp = gpar(fontsize = 12, col = "black", fill = 'black'))
              },
              
              
              # LEGEND
              heatmap_legend_param = list(
                at = c(1,10,20,30),
                #labels = c("Negative",'Positive'),
                title = NULL,
                title_position = 'leftcenter-rot',
                legend_height = unit(2, "cm"),
                adjust = c("right", "top")
              ),
              show_heatmap_legend = TRUE
)

p <- draw(ht,annotation_legend_side = "right",row_dend_side = "left")

topptx(p, filename = 'figures/Figure_Heatmap_C_Comparison_Transcriptome.pptx', 
       width=7, height = 4)

topptx(p, filename = 'figures/Figure_Heatmap_AUC_Comparison_Transcriptome.pptx', 
       width=7, height = 4)

topptx(p, filename = 'figures/Figure_Heatmap_KM_Comparison_Transcriptome.pptx', 
       width=7, height = 4)


## stats

idx <- which(Figure_dataForVolcanoPlot_TCGA$Signature.Count>=3)
idx

sum(Figure_dataForVolcanoPlot_TCGA[idx,]$Significance != 'NS')

genes <- Figure_dataForVolcanoPlot_TCGA$Symbol[idx]
genes


sum(KM_1032_Genes_TCGA$P<0.05, na.rm = T)

idx <- which(KM_1032_Genes_TCGA$Symbol %in% genes)
idx

sum(KM_1032_Genes_TCGA$P[idx]<0.05, na.rm = T)

idx <- which(CoxPH_1032_Genes_TCGA$Symbol %in% genes)
idx

sum(CoxPH_1032_Genes_TCGA$P[idx]<0.05, na.rm = T)


sum(CoxPH_1032_Genes_TCGA$P<0.05, na.rm = T)




