
#########################################################################
###        Evaluation of Survial Analysis Models (For Testing)        ###
#########################################################################

setwd('~/bigdata/PCa/')
source('script/Helper_Functions.R')

# BRB-ArrayTools: https://brb.nci.nih.gov/BRB-ArrayTools/download.html
# randomForestSRC: https://kogalur.github.io/randomForestSRC/theory.html
# superpc: http://statweb.stanford.edu/~tibs/superpc/tutorial.html
# plsRcox: https://fbertran.github.io/plsRcox/

# library(remotes)
# install_github("nyiuab/BhGLM", force=T, build_vignettes=T)


library(Biobase)
library(readxl)
library(survminer)
library(survcomp)

library(survival) # coxph
library(glmnet) # coxnet - lasso, ridge
library(superpc) # pca
library(randomForestSRC) # random forest
library(plsRcox) # pls


library(survivalsvm) # svm
#library(BhGLM)

# CoxPH
# CoxNet-Lasso
# CoxNet-Ridge
# SuperPC
# plsRcox
# RandomForest


# survivalSVM
# Bayesian Losso
# Neural Network Cox
# Dr. Yi


####### TCGA_PRAD #######

datasets <- c('TCGA_PRAD','GSE107299','GSE21034','DKFZ2018','GSE54460','GSE70768','GSE70769',
              'GSE94767','E_MTAB_6128','GSE116918_BCR') # , 'GSE116918_Metastasis'

training.set <- datasets[1]
training.set
eSet <- readRDS(paste0('data/Database/Primary/', training.set, '_eSet.RDS'))
exprData <- exprs(eSet)
phenoData <- pData(eSet)
#View(phenoData)
table(phenoData$sample_type)

keep <- which(phenoData$sample_type=='Primary' | phenoData$sample_type=='Tumor' | phenoData$sample_type=='Tumour')
exprData <- exprData[,keep]
phenoData <- phenoData[keep,]

filter <- which(phenoData$filter=='Duplicate')

if (length(filter)>0) {
  exprData <- exprData[,-filter]
  phenoData <- phenoData[-filter,]
}

filter <- which(duplicated(phenoData$patient_id, incomparables = NA))

if (length(filter)>0) {
  exprData <- exprData[,-filter]
  phenoData <- phenoData[-filter,]
}


filter <- which(is.na(phenoData$bcr_status) | is.na(phenoData$time_to_bcr))

if (length(filter)>0) {
  exprData <- exprData[,-filter]
  phenoData <- phenoData[-filter,]
}


if (training.set %in% c('TCGA_PRAD','GSE54460')) {
  keep <- rowSums(exprData > 0) >= 0.5*ncol(exprData)
  sum(keep)
  exprData <- exprData[keep,]
}


training.pheno <- phenoData
dim(training.pheno)
training.pheno$time_to_bcr[training.pheno$time_to_bcr<=0] <- 0.01

training.geno.all <- as.matrix(t(exprData))
training.geno.all <- scale(training.geno.all)

filter <- which(colSums(is.na(training.geno.all)) !=0)

if (length(filter)>0) {
  training.geno.all <- training.geno.all[,-filter]
}


##############
training.bcr.time <- training.pheno$time_to_bcr
training.bcr.status <-training.pheno$bcr_status

# coxphFun <- function(bcr.time, bcr.status, expr) {
#   
#   coxtest <- coxph(Surv(bcr.time, bcr.status) ~ expr)
#   summcph <- summary(coxtest)
#   
#   coeffs <- c(summcph$coefficients[,2], summcph$conf.int[,3:4], 
#               summcph$coefficients[,5])
#   
#   return (coeffs)
# }
# 
# cox.test <- apply(training.geno.all, 2, function(v) coxphFun(bcr.time, bcr.status, v))
# idx <- which(cox.test[4,] < 0.01)
# cox.genes <- colnames(training.geno.all)[idx]
# cox.genes

cox.list <- list()
for (i in 1:ncol(training.geno.all)) {
  
  print (i)
  expr <- training.geno.all[,i]
  
  if (sd(expr)==0) {
    next
  }
  
  coxtest <- coxph(Surv(training.bcr.time, training.bcr.status) ~ expr)
  summcph <- summary(coxtest)
  
  coeffs <- c(summcph$coefficients[,2], summcph$conf.int[,3:4], 
              summcph$coefficients[,5])
  
  cox.list[[i]] <- coeffs
  
}

cox.table <- do.call(rbind, cox.list)
cox.table <- data.frame(cox.table, stringsAsFactors = F)

colnames(cox.table) <- c('hr','lower95','upper95','p')

cox.table$fdr <- p.adjust(cox.table$p, method = 'BH')

idx <- which(cox.table$p < 0.01)
idx
cox.genes <- colnames(training.geno.all)[idx]
cox.genes
dim(cox.table)

###########################################################################################################

test.sets <- datasets[which(!datasets %in% training.set)]
test.sets

test.set <- test.sets[1]
message (test.set)

eSet <- readRDS(paste0('data/Database/Primary/', test.set, '_eSet.RDS'))
exprData <- exprs(eSet)
phenoData <- pData(eSet)
#View(phenoData)
table(phenoData$sample_type)


keep <- which(phenoData$sample_type=='Primary' | phenoData$sample_type=='Tumor' | phenoData$sample_type=='Tumour')
exprData <- exprData[,keep]
phenoData <- phenoData[keep,]

filter <- which(phenoData$filter=='Duplicate')

if (length(filter)>0) {
  exprData <- exprData[,-filter]
  phenoData <- phenoData[-filter,]
}

filter <- which(duplicated(phenoData$patient_id, incomparables = NA))

if (length(filter)>0) {
  exprData <- exprData[,-filter]
  phenoData <- phenoData[-filter,]
}

filter <- which(is.na(phenoData$bcr_status) | is.na(phenoData$time_to_bcr))

if (length(filter)>0) {
  exprData <- exprData[,-filter]
  phenoData <- phenoData[-filter,]
}


test.pheno <- phenoData
dim(test.pheno)
test.pheno$time_to_bcr[test.pheno$time_to_bcr<=0] <- 0.01

test.geno.all <- as.matrix(t(exprData))
test.geno.all <- scale(test.geno.all)

filter <- which(colSums(is.na(test.geno.all)) !=0)

if (length(filter)>0) {
  test.geno.all <- test.geno.all[,-filter]
}

####################

genes <- Reduce(intersect, list(cox.genes, colnames(training.geno.all), colnames(test.geno.all)))
genes

training.geno <- training.geno.all[,genes]
test.geno <- test.geno.all[,genes]

filter <- which(apply(training.geno, 2, function(v) sum(v==mean(v))==length(v)))
#filter <- which(apply(training.geno, 2, sd)==0)

if (length(filter)>0) {
  training.geno <- training.geno[,-filter]
  test.geno <- test.geno[,-filter]
}

# 
# 
# 
# model <- 'CoxPH'
# print (model)
# 
# training.bcr.time <- training.pheno$time_to_bcr
# training.bcr.status <-training.pheno$bcr_status
# 
# test.bcr.time <- test.pheno$time_to_bcr
# test.bcr.status <-  test.pheno$bcr_status
# 
# training.surv.data <- data.frame(training.geno, bcr.time=training.bcr.time, bcr.status=training.bcr.status)
# 
# multi.var <- paste0(genes, collapse = '+')
# 
# fml <- as.formula(paste0('Surv(bcr.time, bcr.status) ~ ', multi.var))
# 
# coxtest <- coxph(fml, data = training.surv.data)
# summcph <- summary(coxtest)
# coeffs <- summcph$coefficients[,1]
# 
# training.risk.score <- predict(coxtest, data.frame(training.geno), type="lp", se.fit=FALSE) # type='risk'
# training.risk.score
# 
# training.risk.threshold <- median(training.risk.score, na.rm = T)
# training.risk.threshold
# 
# test.risk.score <- predict(coxtest, data.frame(test.geno), type="lp", se.fit=FALSE)
# test.risk.score
# 
# test.risk.threshold <- median(test.risk.score, na.rm = T)
# test.risk.threshold
# 
# training.risk.group <- test.risk.score > training.risk.threshold
# training.risk.group
# 
# test.risk.group <- test.risk.score > test.risk.threshold
# test.risk.group
# 
# riskScore <- test.risk.score
# trainingRiskGroup <- training.risk.group
# testRiskGroup <- test.risk.group
# 
# riskThreshold <- median(riskScore, na.rm=T)
# riskGroup <- riskScore > riskThreshold
# 
# res <- data.frame(risk.score=riskScore, risk.group=riskGroup, training.risk.group=trainingRiskGroup,
#                   test.risk.group=testRiskGroup, bcr.time=test.pheno$time_to_bcr, bcr.status=test.pheno$bcr_status)
# res

############################## CoxLasso

alpha <- 1
model <- 'CoxLasso'

training.bcr.time <- training.pheno$time_to_bcr
training.bcr.status <-training.pheno$bcr_status

test.bcr.time <- test.pheno$time_to_bcr
test.bcr.status <-  test.pheno$bcr_status

training.surv.data <- data.frame(training.geno, bcr.time=training.bcr.time, bcr.status=training.bcr.status)

set.seed(777)
cv.fit <- cv.glmnet(training.geno, Surv(training.bcr.time, training.bcr.status), family="cox", maxit = 1000,
                    alpha=alpha)

coeffs <- coef(cv.fit, s = cv.fit$lambda.min)
coeffs <- as.numeric(coeffs)

message(paste0(sum(coeffs!=0),'/',length(genes)))

training.risk.score <- predict(cv.fit, s=cv.fit$lambda.min, newx = training.geno, type="link") # response: reltive risk; link: linear prediction
training.risk.score
training.risk.score <- training.risk.score[,1]

training.risk.threshold <- median(training.risk.score, na.rm = T)
training.risk.threshold

test.risk.score <- predict(cv.fit, s=cv.fit$lambda.min, newx = test.geno, type="link") # response: reltive risk; link: linear prediction
test.risk.score
test.risk.score <- test.risk.score[,1]

test.risk.threshold <- median(test.risk.score, na.rm = T)
test.risk.threshold

training.risk.group <- test.risk.score > training.risk.threshold
training.risk.group

test.risk.group <- test.risk.score > test.risk.threshold
test.risk.group

riskScore <- test.risk.score
trainingRiskGroup <- training.risk.group
testRiskGroup <- test.risk.group

riskThreshold <- median(riskScore, na.rm=T)
riskGroup <- riskScore > riskThreshold

res <- data.frame(risk.score=riskScore, risk.group=riskGroup, training.risk.group=trainingRiskGroup,
                  test.risk.group=testRiskGroup, bcr.time=test.pheno$time_to_bcr, bcr.status=test.pheno$bcr_status)
res

#write.csv(res, file=paste0('report/Survival/TrainTestScale/Transcriptome/Survival_',
#                           model, '_', signature.name, '_', training.set, '_', test.set, '_TrainTest_Transcriptome.csv'),
#          quote = F)



############################## CoxRidge

alpha <- 0
model <- 'CoxRidge'

training.bcr.time <- training.pheno$time_to_bcr
training.bcr.status <-training.pheno$bcr_status

test.bcr.time <- test.pheno$time_to_bcr
test.bcr.status <-  test.pheno$bcr_status

training.surv.data <- data.frame(training.geno, bcr.time=training.bcr.time, bcr.status=training.bcr.status)

set.seed(777)
cv.fit <- cv.glmnet(training.geno, Surv(training.bcr.time, training.bcr.status), family="cox", maxit = 1000,
                    alpha=alpha)

coeffs <- coef(cv.fit, s = cv.fit$lambda.min)
coeffs <- as.numeric(coeffs)

message(paste0(sum(coeffs!=0),'/',length(genes)))

training.risk.score <- predict(cv.fit, s=cv.fit$lambda.min, newx = training.geno, type="link") # response: reltive risk; link: linear prediction
training.risk.score
training.risk.score <- training.risk.score[,1]

training.risk.threshold <- median(training.risk.score, na.rm = T)
training.risk.threshold

test.risk.score <- predict(cv.fit, s=cv.fit$lambda.min, newx = test.geno, type="link") # response: reltive risk; link: linear prediction
test.risk.score
test.risk.score <- test.risk.score[,1]

test.risk.threshold <- median(test.risk.score, na.rm = T)
test.risk.threshold

training.risk.group <- test.risk.score > training.risk.threshold
training.risk.group

test.risk.group <- test.risk.score > test.risk.threshold
test.risk.group

riskScore <- test.risk.score
trainingRiskGroup <- training.risk.group
testRiskGroup <- test.risk.group

riskThreshold <- median(riskScore, na.rm=T)
riskGroup <- riskScore > riskThreshold

res <- data.frame(risk.score=riskScore, risk.group=riskGroup, training.risk.group=trainingRiskGroup,
                  test.risk.group=testRiskGroup, bcr.time=test.pheno$time_to_bcr, bcr.status=test.pheno$bcr_status)
res



##################################### SuperPC

# threshold <- 0 # 0.1, 0.3
# model <- paste0('SuperPC',threshold)

threshold <- 1
model <- 'SuperPC'

n.components <- 1

training.bcr.time <- training.pheno$time_to_bcr
training.bcr.status <-training.pheno$bcr_status

test.bcr.time <- test.pheno$time_to_bcr
test.bcr.status <-  test.pheno$bcr_status

#training.surv.data <- data.frame(training.geno, bcr.time=training.bcr.time, bcr.status=training.bcr.status)

training.data.pca <-list(x=t(training.geno), y=training.bcr.time, censoring.status= training.bcr.status, featurenames=genes)
test.data.pca <- list(x=t(test.geno), y=test.bcr.time, censoring.status= test.bcr.status, featurenames=genes)

pca.fit <- superpc.train(training.data.pca, type="survival")

message(paste0(sum(abs(pca.fit$feature.scores)>threshold),'/',length(genes)))

training.risk.score <- superpc.predict(pca.fit, training.data.pca, training.data.pca,
                                       threshold=threshold, n.components=n.components, prediction.type="continuous")

training.risk.score <- training.risk.score$v.pred[,1]
names(training.risk.score) <- rownames(training.geno)

training.risk.threshold <- median(training.risk.score, na.rm = T)

test.risk.score <- superpc.predict(pca.fit, training.data.pca, test.data.pca,
                                   threshold=threshold, n.components=n.components, prediction.type="continuous")

test.risk.score <- test.risk.score$v.pred[,1]
names(test.risk.score) <- rownames(test.geno)

test.risk.threshold <- median(test.risk.score, na.rm = T)
test.risk.threshold

training.risk.group <- test.risk.score > training.risk.threshold
training.risk.group

test.risk.group <- test.risk.score > test.risk.threshold
test.risk.group

riskScore <- test.risk.score
trainingRiskGroup <- training.risk.group
testRiskGroup <- test.risk.group

riskThreshold <- median(riskScore, na.rm=T)
riskGroup <- riskScore > riskThreshold

res <- data.frame(risk.score=riskScore, risk.group=riskGroup, training.risk.group=trainingRiskGroup,
                  test.risk.group=testRiskGroup, bcr.time=test.pheno$time_to_bcr, bcr.status=test.pheno$bcr_status)
res



##################################### plsRcox

ncomps <- 2
model <- 'plsRcox'


training.bcr.time <- training.pheno$time_to_bcr
training.bcr.status <-training.pheno$bcr_status

test.bcr.time <- test.pheno$time_to_bcr
test.bcr.status <-  test.pheno$bcr_status

training.surv.data <- data.frame(training.geno, bcr.time=training.bcr.time, bcr.status=training.bcr.status)


pls.fit <- plsRcox(training.geno,time=training.bcr.time,event=training.bcr.status, nt=ncomps)

training.risk.score <- predict(pls.fit, newdata = training.geno, type="lp", comps= ncomps, se.fit=FALSE) # type='risk'
training.risk.score

names(training.risk.score) <- rownames(training.geno)

training.risk.threshold <- median(training.risk.score, na.rm = T)
training.risk.threshold

test.risk.score <- predict(pls.fit, newdata = test.geno, type="lp", comps= ncomps, se.fit=FALSE)
test.risk.score
#risk.threshold <- median(test.risk.score, na.rm = T)
#risk.threshold

names(test.risk.score) <- rownames(test.geno)

test.risk.threshold <- median(test.risk.score, na.rm = T)
test.risk.threshold

training.risk.group <- test.risk.score > training.risk.threshold
training.risk.group

test.risk.group <- test.risk.score > test.risk.threshold
test.risk.group

riskScore <- test.risk.score
trainingRiskGroup <- training.risk.group
testRiskGroup <- test.risk.group

riskThreshold <- median(riskScore, na.rm=T)
riskGroup <- riskScore > riskThreshold

res <- data.frame(risk.score=riskScore, risk.group=riskGroup, training.risk.group=trainingRiskGroup,
                  test.risk.group=testRiskGroup, bcr.time=test.pheno$time_to_bcr, bcr.status=test.pheno$bcr_status)
res



##################################### RandomForest

ntree <- 100
model <- 'RandomForest'


training.bcr.time <- training.pheno$time_to_bcr
training.bcr.status <-training.pheno$bcr_status

test.bcr.time <- test.pheno$time_to_bcr
test.bcr.status <-  test.pheno$bcr_status

training.surv.data <- data.frame(training.geno, bcr.time=training.bcr.time, bcr.status=training.bcr.status)

set.seed(777)
rf.fit <- rfsrc(Surv(bcr.time, bcr.status) ~ ., data = training.surv.data, ntree = ntree)

training.risk.score <- predict(rf.fit, data.frame(training.geno))
training.risk.score <- training.risk.score$predicted
names(training.risk.score) <- rownames(training.geno)

training.risk.threshold <- median(training.risk.score, na.rm = T)
training.risk.threshold

test.risk.score <- predict(rf.fit, data.frame(test.geno)) ###
test.risk.score <- test.risk.score$predicted
names(test.risk.score) <- rownames(test.geno)

test.risk.threshold <- median(test.risk.score, na.rm = T)
test.risk.threshold

training.risk.group <- test.risk.score > training.risk.threshold
training.risk.group

test.risk.group <- test.risk.score > test.risk.threshold
test.risk.group

riskScore <- test.risk.score
trainingRiskGroup <- training.risk.group
testRiskGroup <- test.risk.group

riskThreshold <- median(riskScore, na.rm=T)
riskGroup <- riskScore > riskThreshold

res <- data.frame(risk.score=riskScore, risk.group=riskGroup, training.risk.group=trainingRiskGroup,
                  test.risk.group=testRiskGroup, bcr.time=test.pheno$time_to_bcr, bcr.status=test.pheno$bcr_status)
res


#####################################################################################################

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

