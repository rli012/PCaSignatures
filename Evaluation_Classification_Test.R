
#########################################################################
###         Evaluation of Classification Models (For Testing)         ###
#########################################################################

setwd('~/bigdata/PCa/')
source('script/Helper_Functions.R')

library(Biobase)
library(readxl)
library(caret)
library(glmnet)
library(kernlab)
library(pls)
library(randomForest)
library(MASS) # lda
library(xgboost)
library(pROC)
library(ROCR)

options(stringsAsFactors = F)


####### TCGA_PRAD #######

datasets <- c('TCGA_PRAD','GSE107299','GSE21034','DKFZ2018','GSE54460','GSE70768','GSE70769',
              'GSE94767','E_MTAB_6128','GSE116918_BCR') # , 'GSE116918_Metastasis'

dataset <- datasets[10]

eSet <- readRDS(paste0('data/Database/Primary/', dataset, '_eSet.RDS'))
exprData <- exprs(eSet)
phenoData <- pData(eSet)
#View(phenoData)
table(phenoData$sample_type)

keep <- which(phenoData$sample_type=='Primary' | phenoData$sample_type=='Tumor' | phenoData$sample_type=='Tumour')
exprData <- exprData[,keep]
phenoData <- phenoData[keep,]

filter <- which(duplicated(phenoData$patient_id, incomparables = NA))

if (length(filter)>0) {
  exprData <- exprData[,-filter]
  phenoData <- phenoData[-filter,]
}

#phenoData$bcr_status_3yr <- getEventFun(n=3, time.to.event=phenoData$time_to_bcr, status=phenoData$bcr_status)
phenoData$risk_group <- getEventFun(n=5, time.to.event=phenoData$time_to_bcr, status=phenoData$bcr_status)
#phenoData$bcr_status_10yr <- getEventFun(n=10, time.to.event=phenoData$time_to_bcr, status=phenoData$bcr_status)

#table(phenoData$bcr_status_3yr)
table(phenoData$risk_group)
#table(phenoData$bcr_status_10yr)

filter <- which(is.na(phenoData$risk_group))

if (length(filter)>0) {
  exprData <- exprData[,-filter]
  phenoData <- phenoData[-filter,]
}

pheno <- as.matrix(phenoData$risk_group, drop=FALSE) ### 5 Year
geno <- as.matrix(t(exprData))
geno <- scale(geno)

set.seed(777)
nfold <- 10
foldid <- generateCVFold(n=nrow(geno), nfold=nfold)
#saveRDS(foldid, file=paste0('data/foldid/Foldid_', dataset, '_bcr_5yr.RDS'))


# Published signatures
signatures <- c('Agell','Bibikova','Bismar','Decipher','Ding','Glinsky','Irshad',
                'Jennifer','Jia','Kamoun','Li','Long','Luca','Mo','Nakagawa','Olmos',
                'Oncotype','Penney','Planche','Prolaris','Ramaswamy','Ramos_Montoya',
                'Ross_Adams','Ross_Robert','Sharma','Talantov','Varambally','Wu','Yang',
                'Yu')

signature.name <- signatures[11]

signature <- read_xlsx(path = 'data/PCa_Prognosis_Signatures.xlsx', sheet=signature.name)
signature.genes <- signature$Ensembl

# Random signatures
#random.signatures <- readRDS(file='data/Random_Signatures.RDS')
#signature.genes <- random.signatures[[1]]


models <- c('glmnet','svmLinear','svmRadial','svmPoly','rf','pls','lda','xgbLinear','xgbTree') # dnn
model <- models[1]

res <- runCaretCV(geno = geno, pheno = pheno, foldid = foldid, signature = signature.genes, model = model)
res


samples <- rownames(res)
res$risk.score <- res$yprob
res$bcr.time <- phenoData[samples,]$time_to_bcr
res$bcr.status <- phenoData[samples,]$bcr_status

risk.threshold <- median(res$risk.score, na.rm = T)
res$risk.group <- res$risk.score > risk.threshold

###
risk.score <- res$risk.score
risk.group <- res$risk.group
bcr.time <- res$bcr.time
bcr.status <- res$bcr.status

############ AUC for Classification

table(pred=res$ypred, true=res$yobs)

a <- sum(res$yobs=='E1' & res$ypred=='E1')
b <- sum(res$yobs=='E0' & res$ypred=='E1')
c <- sum(res$yobs=='E1' & res$ypred=='E0')
d <- sum(res$yobs=='E0' & res$ypred=='E0')

sensitivity <- a/(a+c)
specificity <- d/(b+d)

accuracy <- (a+d)/(a+b+c+d)

sensitivity
specificity
accuracy


roc.test <- roc(res$yobs, res$yprob, plot=TRUE, ci=TRUE, auc=TRUE)
ci.auc <- roc.test$ci
ci.auc
auc <- ci.auc[2]
auc.ci.lower95 <- ci.auc[1]
auc.ci.upper95 <- ci.auc[3]

auc <- format(auc, digits = 2, nsmall=2)
auc.ci.lower95 <- format(auc.ci.lower95, digits = 2, nsmall=2)
auc.ci.upper95 <- format(auc.ci.upper95, digits = 2, nsmall=2)

#FPR <- 1-roc.test$specificities
#TPR <- roc.test$sensitivities

pred <- prediction(res$yprob, res$yobs)
pred
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
perf

FPR <- perf@x.values[[1]]
TPR <- perf@y.values[[1]]
#FPR

df <- data.frame(FPR,TPR)

auc.test <- wilcox.test(FPR, TPR, alternative = 'two.sided')
auc.test$p.value
pvalue <- formatC(auc.test$p.value, format = 'e', digits = 2)
pvalue

ggplot(df,aes(x=FPR,y=TPR))+geom_line(size = 1, alpha = 1,color='red')+
  labs(x = "1-Specificity",y = "Sensitivity")+ 
  #scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+
  geom_abline(intercept = 0, slope = 1) +
  #geom_segment(x=0,y=0,xend=1,yend=1, color='darkgreen') + xlim(0,1) + ylim(0,1) +
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black'),
        panel.background = element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) +
  theme(strip.text.x = element_text(size = 12, colour = "black", angle=0)) +
  ggplot2::annotate("text", 
                    x = 0.6, y = 0.125, # x and y coordinates of the text
                    label = paste('AUC=',auc, ' (95% CI: ',auc.ci.lower95, '-', auc.ci.upper95, ')',sep=''), size = 5) +
  ggplot2::annotate("text", 
                    x = 0.6, y = 0.25, # x and y coordinates of the text
                    label = paste0('P=',pvalue), size = 5)


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

td.auc <- roc$AUC
td.auc

############ CoxPH

coxtest <- coxph(Surv(bcr.time, bcr.status) ~ risk.score)
summcph <- summary(coxtest)

coeffs <- c(summcph$coefficients[,2], summcph$conf.int[,3:4], 
            summcph$coefficients[,5])

coeffs


############## KM Plot

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
