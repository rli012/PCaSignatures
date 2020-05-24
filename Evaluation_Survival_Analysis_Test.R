
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
library(BhGLM)

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

dataset <- datasets[1]

eSet <- readRDS(paste0('data/Database/Primary/', dataset, '_eSet.RDS'))
exprData <- exprs(eSet)
phenoData <- pData(eSet)
#View(phenoData)
table(phenoData$sample_type)
keep <- which(phenoData$sample_type=='Primary')
keep

exprData <- exprData[,keep]
phenoData <- phenoData[keep,]

filter <- which(is.na(phenoData$bcr_status) | is.na(phenoData$time_to_bcr))
filter
phenoData <- phenoData[-filter,]
exprData <- exprData[,-filter]

pheno <- phenoData
pheno$time_to_bcr[pheno$time_to_bcr==0] <- 0.01

geno <- as.matrix(t(exprData))
geno <- scale(geno)

samples <- rownames(geno)
samples

set.seed(777)
nfold <- 10
foldid <- generateCVFold(nsam=nrow(geno), nfold=nfold)
#saveRDS(foldid, file=paste0('data/foldid/Foldid_', dataset, '_bcr_5yr.RDS'))
foldid


# Published signatures
signatures <- c('Agell','Bibikova','Bismar','Decipher','Ding','Glinsky','Irshad',
                'Jennifer','Jia','Kamoun','Li','Long','Luca','Mo','Nakagawa','Olmos',
                'Oncotype','Penney','Planche','Prolaris','Ramaswamy','Ramos_Montoya',
                'Ross_Adams','Ross_Robert','Sharma','Talantov','Varambally','Wu','Yang',
                'Yu')


### Test a single signature
signature.name <- signatures[1]

signature <- read_xlsx(path = 'data/PCa_Prognosis_Signatures.xlsx', sheet=signature.name)
signature <- signature$Ensembl
signature

signature <- signature[which(signature!='NA')]
signature

genes <- signature

#genes <- sample(colnames(geno),10)
#genes


############################## CoxPH

riskScore <- c()
trainingRiskGroup <- c()
testRiskGroup <- c()


cv <- list()

for (k in 1:nfold) {
  
  training.idx <- which(foldid!=k)
  test.idx <- which(foldid==k)
  
  training.geno <- geno[training.idx,genes]
  test.geno <- geno[test.idx,genes]
  
  filter <- which(apply(training.geno, 2, function(v) sum(v==mean(v))==length(v)))
  #filter <- which(apply(training.geno, 2, sd)==0)
  
  if (length(filter)>0) {
    training.geno <- training.geno[,-filter]
    test.geno <- test.geno[,-filter]
  }
  

  training.bcr.time <- pheno$time_to_bcr[training.idx]
  training.bcr.status <-  pheno$bcr_status[training.idx]
  
  test.bcr.time <- pheno$time_to_bcr[test.idx]
  test.bcr.status <-  pheno$bcr_status[test.idx]
  
  training.surv.data <- data.frame(training.geno, bcr.time=training.bcr.time, bcr.status=training.bcr.status)
  
  multi.var <- paste0(genes, collapse = '+')
  
  fml <- as.formula(paste0('Surv(bcr.time, bcr.status) ~ ', multi.var))
  
  coxtest <- coxph(fml, data = training.surv.data)
  summcph <- summary(coxtest)
  coeffs <- summcph$coefficients[,1]
  
  # risk.score1 <- coxtest$linear.predictors # NA removed
  
  #risk.score2 <- colSums(apply(training.geno, 1, function(v) v*coeffs))
  #risk.score2
  
  #risk.score3 <- predict(coxtest, data.frame(training.geno), type="lp", se.fit=FALSE)
  #risk.score3
  
  #risk.score2-risk.score3 # not exactly the same
  
  #risk.score <- predict(coxtest, data.frame(training.geno), type="risk", se.fit=FALSE)
  #risk.score
  
  training.risk.score <- predict(coxtest, data.frame(training.geno), type="risk", se.fit=FALSE)
  training.risk.score
  
  training.risk.threshold <- median(training.risk.score, na.rm = T)
  training.risk.threshold
  
  test.risk.score <- predict(coxtest, data.frame(test.geno), type="risk", se.fit=FALSE)
  test.risk.score
  #risk.threshold <- median(test.risk.score, na.rm = T)
  #risk.threshold
  
  test.risk.threshold <- median(test.risk.score, na.rm = T)
  test.risk.threshold
  
  
  training.risk.group <- test.risk.score > training.risk.threshold
  training.risk.group
  
  test.risk.group <- test.risk.score > test.risk.threshold
  test.risk.group
  
  tr <- paste0('[Train] ', training.risk.score)
  names(tr) <- names(training.risk.score)
  
  test <- paste0('[Test] ', test.risk.score)
  names(test) <- names(test.risk.score)
  
  cv[[k]] <- c(tr, test)[samples]
  
  riskScore <- c(riskScore, test.risk.score)
  
  trainingRiskGroup <- c(trainingRiskGroup, training.risk.group)
  testRiskGroup <- c(testRiskGroup, test.risk.group)
  
}

cv <- do.call(cbind, cv)
colnames(cv) <- paste0('CV',1:10)

riskScore <- riskScore[samples]
trainingRiskGroup <- trainingRiskGroup[samples]
testRiskGroup <- testRiskGroup[samples]

riskThreshold <- median(riskScore, na.rm=T)
riskGroup <- riskScore > riskThreshold

res <- data.frame(risk.score=riskScore, risk.group=riskGroup, training.risk.group=trainingRiskGroup, 
                  test.risk.group=testRiskGroup, bcr.time=pheno$time_to_bcr, bcr.status=pheno$bcr_status, cv)


############################## Elastic Net (CoxNet) - LASSO

riskScore <- c()
trainingRiskGroup <- c()
testRiskGroup <- c()

cv <- list()

for (k in 1:nfold) {
  
  training.idx <- which(foldid!=k)
  test.idx <- which(foldid==k)
  
  training.geno <- geno[training.idx,genes]
  test.geno <- geno[test.idx,genes]
  
  filter <- which(apply(training.geno, 2, function(v) sum(v==mean(v))==length(v)))
  #filter <- which(apply(training.geno, 2, sd)==0)
  
  if (length(filter)>0) {
    training.geno <- training.geno[,-filter]
    test.geno <- test.geno[,-filter]
  }
  
  training.bcr.time <- pheno$time_to_bcr[training.idx]
  training.bcr.status <-  pheno$bcr_status[training.idx]
  
  test.bcr.time <- pheno$time_to_bcr[test.idx]
  test.bcr.status <-  pheno$bcr_status[test.idx]
  
  training.surv.data <- data.frame(training.geno, bcr.time=training.bcr.time, bcr.status=training.bcr.status)
  
  set.seed(777)
  cv.fit <- cv.glmnet(training.geno, Surv(training.bcr.time, training.bcr.status), family="cox", alpha = 1, maxit = 1000)
  
  coeffs <- coef(cv.fit, s = cv.fit$lambda.min)
  coeffs <- as.numeric(coeffs)
  
  training.risk.score <- predict(cv.fit, s=cv.fit$lambda.min, newx = training.geno, type="response") # response: reltive risk; link: linear prediction
  training.risk.score
  training.risk.score <- training.risk.score[,1]
  
  training.risk.threshold <- median(training.risk.score, na.rm = T)
  training.risk.threshold
  
  test.risk.score <- predict(cv.fit, s=cv.fit$lambda.min, newx = test.geno, type="response") # response: reltive risk; link: linear prediction
  test.risk.score
  test.risk.score <- test.risk.score[,1]
  
  test.risk.threshold <- median(test.risk.score, na.rm = T)
  test.risk.threshold
  
  
  training.risk.group <- test.risk.score > training.risk.threshold
  training.risk.group
  
  test.risk.group <- test.risk.score > test.risk.threshold
  test.risk.group
  
  tr <- paste0('[Train] ', training.risk.score)
  names(tr) <- names(training.risk.score)
  
  test <- paste0('[Test] ', test.risk.score)
  names(test) <- names(test.risk.score)
  
  cv[[k]] <- c(tr, test)[samples]
  
  riskScore <- c(riskScore, test.risk.score)
  
  trainingRiskGroup <- c(trainingRiskGroup, training.risk.group)
  testRiskGroup <- c(testRiskGroup, test.risk.group)
  
}

cv <- do.call(cbind, cv)
colnames(cv) <- paste0('CV',1:10)
cv

riskScore <- riskScore[samples]
trainingRiskGroup <- trainingRiskGroup[samples]
testRiskGroup <- testRiskGroup[samples]

riskThreshold <- median(riskScore, na.rm=T)
riskGroup <- riskScore > riskThreshold

res <- data.frame(risk.score=riskScore, risk.group=riskGroup, training.risk.group=trainingRiskGroup, 
                  test.risk.group=testRiskGroup, bcr.time=pheno$time_to_bcr, bcr.status=pheno$bcr_status, cv)


############################## Elastic Net (CoxNet) - Ridge


riskScore <- c()
trainingRiskGroup <- c()
testRiskGroup <- c()

cv <- list()

for (k in 1:nfold) {
  
  training.idx <- which(foldid!=k)
  test.idx <- which(foldid==k)
  
  training.geno <- geno[training.idx,genes]
  test.geno <- geno[test.idx,genes]
  
  filter <- which(apply(training.geno, 2, function(v) sum(v==mean(v))==length(v)))
  #filter <- which(apply(training.geno, 2, sd)==0)
  
  if (length(filter)>0) {
    training.geno <- training.geno[,-filter]
    test.geno <- test.geno[,-filter]
  }
  
  training.bcr.time <- pheno$time_to_bcr[training.idx]
  training.bcr.status <-  pheno$bcr_status[training.idx]
  
  test.bcr.time <- pheno$time_to_bcr[test.idx]
  test.bcr.status <-  pheno$bcr_status[test.idx]
  
  training.surv.data <- data.frame(training.geno, bcr.time=training.bcr.time, bcr.status=training.bcr.status)
  
  set.seed(777)
  cv.fit <- cv.glmnet(training.geno, Surv(training.bcr.time, training.bcr.status), family="cox", alpha = 0, maxit = 1000)
  
  coeffs <- coef(cv.fit, s = cv.fit$lambda.min)
  coeffs <- as.numeric(coeffs)
  
  training.risk.score <- predict(cv.fit, s=cv.fit$lambda.min, newx = training.geno, type="response") # response: reltive risk; link: linear prediction
  training.risk.score
  training.risk.score <- training.risk.score[,1]
  
  training.risk.threshold <- median(training.risk.score, na.rm = T)
  training.risk.threshold
  
  test.risk.score <- predict(cv.fit, s=cv.fit$lambda.min, newx = test.geno, type="response") # response: reltive risk; link: linear prediction
  test.risk.score
  test.risk.score <- test.risk.score[,1]
  
  test.risk.threshold <- median(test.risk.score, na.rm = T)
  test.risk.threshold
  
  
  training.risk.group <- test.risk.score > training.risk.threshold
  training.risk.group
  
  test.risk.group <- test.risk.score > test.risk.threshold
  test.risk.group
  
  tr <- paste0('[Train] ', training.risk.score)
  names(tr) <- names(training.risk.score)
  
  test <- paste0('[Test] ', test.risk.score)
  names(test) <- names(test.risk.score)
  
  cv[[k]] <- c(tr, test)[samples]
  
  riskScore <- c(riskScore, test.risk.score)
  
  trainingRiskGroup <- c(trainingRiskGroup, training.risk.group)
  testRiskGroup <- c(testRiskGroup, test.risk.group)
  
}

cv <- do.call(cbind, cv)
colnames(cv) <- paste0('CV',1:10)
cv

riskScore <- riskScore[samples]
trainingRiskGroup <- trainingRiskGroup[samples]
testRiskGroup <- testRiskGroup[samples]

riskThreshold <- median(riskScore, na.rm=T)
riskGroup <- riskScore > riskThreshold

res <- data.frame(risk.score=riskScore, risk.group=riskGroup, training.risk.group=trainingRiskGroup, 
                  test.risk.group=testRiskGroup, bcr.time=pheno$time_to_bcr, bcr.status=pheno$bcr_status, cv)


##################################### superpc

threshold <- 0.3

riskScore <- c()
trainingRiskGroup <- c()
testRiskGroup <- c()

cv <- list()

for (k in 1:nfold) {
  
  training.idx <- which(foldid!=k)
  test.idx <- which(foldid==k)
  
  training.geno <- geno[training.idx,genes]
  test.geno <- geno[test.idx,genes]
  
  filter <- which(apply(training.geno, 2, function(v) sum(v==mean(v))==length(v)))
  #filter <- which(apply(training.geno, 2, sd)==0)
  
  if (length(filter)>0) {
    training.geno <- training.geno[,-filter]
    test.geno <- test.geno[,-filter]
  }
  
  training.bcr.time <- pheno$time_to_bcr[training.idx]
  training.bcr.status <-  pheno$bcr_status[training.idx]
  
  test.bcr.time <- pheno$time_to_bcr[test.idx]
  test.bcr.status <-  pheno$bcr_status[test.idx]
  
  training.data.pca <-list(x=t(training.geno), y=training.bcr.time, censoring.status= training.bcr.status, featurenames=genes)
  test.data.pca <- list(x=t(test.geno), y=test.bcr.time, censoring.status= test.bcr.status, featurenames=genes)
  
  set.seed(777)
  pca.fit <- superpc.train(training.data.pca, type="survival")
  
  training.risk.score <- superpc.predict(pca.fit, training.data.pca, training.data.pca, 
                                         threshold=threshold, n.components=1, prediction.type="continuous")
  
  training.risk.score <- training.risk.score$v.pred[,1]
  names(training.risk.score) <- rownames(training.geno)
  
  test.risk.score <- superpc.predict(pca.fit, training.data.pca, test.data.pca, 
                                     threshold=threshold, n.components=1, prediction.type="continuous")
  
  test.risk.score <- test.risk.score$v.pred[,1]
  names(test.risk.score) <- rownames(test.geno)
  
  #training.surv.data <- data.frame(training.geno, bcr.time=training.bcr.time, bcr.status=training.bcr.status)
  
  test.risk.threshold <- median(test.risk.score, na.rm = T)
  test.risk.threshold
  
  
  training.risk.group <- test.risk.score > training.risk.threshold
  training.risk.group
  
  test.risk.group <- test.risk.score > test.risk.threshold
  test.risk.group
  
  tr <- paste0('[Train] ', training.risk.score)
  names(tr) <- names(training.risk.score)
  
  test <- paste0('[Test] ', test.risk.score)
  names(test) <- names(test.risk.score)
  
  cv[[k]] <- c(tr, test)[samples]
  
  riskScore <- c(riskScore, test.risk.score)
  
  trainingRiskGroup <- c(trainingRiskGroup, training.risk.group)
  testRiskGroup <- c(testRiskGroup, test.risk.group)
  
}

cv <- do.call(cbind, cv)
colnames(cv) <- paste0('CV',1:10)
cv

riskScore <- riskScore[samples]
trainingRiskGroup <- trainingRiskGroup[samples]
testRiskGroup <- testRiskGroup[samples]

riskThreshold <- median(riskScore, na.rm=T)
riskGroup <- riskScore > riskThreshold

res <- data.frame(risk.score=riskScore, risk.group=riskGroup, training.risk.group=trainingRiskGroup, 
                  test.risk.group=testRiskGroup, bcr.time=pheno$time_to_bcr, bcr.status=pheno$bcr_status, cv)


##################################### plsRcox

ncomps <- 2

riskScore <- c()
trainingRiskGroup <- c()
testRiskGroup <- c()

cv <- list()

for (k in 1:nfold) {
  
  training.idx <- which(foldid!=k)
  test.idx <- which(foldid==k)
  
  training.geno <- geno[training.idx,genes]
  test.geno <- geno[test.idx,genes]
  
  filter <- which(apply(training.geno, 2, function(v) sum(v==mean(v))==length(v)))
  #filter <- which(apply(training.geno, 2, sd)==0)
  
  if (length(filter)>0) {
    training.geno <- training.geno[,-filter]
    test.geno <- test.geno[,-filter]
  }
  
  training.bcr.time <- pheno$time_to_bcr[training.idx]
  training.bcr.status <-  pheno$bcr_status[training.idx]
  
  test.bcr.time <- pheno$time_to_bcr[test.idx]
  test.bcr.status <-  pheno$bcr_status[test.idx]
  
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
  
  tr <- paste0('[Train] ', training.risk.score)
  names(tr) <- names(training.risk.score)
  
  test <- paste0('[Test] ', test.risk.score)
  names(test) <- names(test.risk.score)
  
  cv[[k]] <- c(tr, test)[samples]
  
  riskScore <- c(riskScore, test.risk.score)
  
  trainingRiskGroup <- c(trainingRiskGroup, training.risk.group)
  testRiskGroup <- c(testRiskGroup, test.risk.group)
  
}

cv <- do.call(cbind, cv)
cv
colnames(cv) <- paste0('CV',1:10)

riskScore <- riskScore[samples]
trainingRiskGroup <- trainingRiskGroup[samples]
testRiskGroup <- testRiskGroup[samples]

riskThreshold <- median(riskScore, na.rm=T)
riskGroup <- riskScore > riskThreshold

res <- data.frame(risk.score=riskScore, risk.group=riskGroup, training.risk.group=trainingRiskGroup,
                  test.risk.group=testRiskGroup, bcr.time=pheno$time_to_bcr, bcr.status=pheno$bcr_status, cv)


# coxsplsDR and cv.coxsplsDR
#(cox_splsDR_fit=coxsplsDR(X_train_micro,Y_train_micro,C_train_micro,ncomp=6,eta=.5))
#(cv.coxsplsDR.res=cv.coxsplsDR(list(x=X_train_micro,time=Y_train_micro, status=C_train_micro),nt=10,eta=.1))


# coxDKsplsDR and cv.coxDKsplsDR
#(cox_DKsplsDR_fit=coxDKsplsDR(X_train_micro,Y_train_micro,C_train_micro,ncomp=6, validation="CV",eta=.5))
#(cv.coxDKsplsDR.res=cv.coxDKsplsDR(list(x=X_train_micro,time=Y_train_micro, status=C_train_micro),nt=10,eta=.1))


# coxDKplsDR and cv.coxDKplsDR

# coxpls and cv.coxpls

# coxplsDR and cv.coxplsDR

# DKplsRcox
#DKplsRcox(X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5)

# larsDR and cv.larsDR


# plsRcox and cv.plsRcox
#plsRcox(Xplan=~X_train_micro,time=Y_train_micro,event=C_train_micro,nt=5,sparse=TRUE, alpha.pvals.expli=.15)
#(cv.plsRcox.res=cv.plsRcox(list(x=X_train_micro,time=Y_train_micro,status=C_train_micro),nt=10,verbose = FALSE))


# autoplsRcox and cv.autoplsRcox


##################################### randomForestSRC

ntree <- 100

riskScore <- c()
trainingRiskGroup <- c()
testRiskGroup <- c()

cv <- list()

for (k in 1:nfold) {
  
  training.idx <- which(foldid!=k)
  test.idx <- which(foldid==k)
  
  training.geno <- geno[training.idx,genes]
  test.geno <- geno[test.idx,genes]
  
  filter <- which(apply(training.geno, 2, function(v) sum(v==mean(v))==length(v)))
  #filter <- which(apply(training.geno, 2, sd)==0)
  
  if (length(filter)>0) {
    training.geno <- training.geno[,-filter]
    test.geno <- test.geno[,-filter]
  }
  
  training.bcr.time <- pheno$time_to_bcr[training.idx]
  training.bcr.status <-  pheno$bcr_status[training.idx]
  
  test.bcr.time <- pheno$time_to_bcr[test.idx]
  test.bcr.status <-  pheno$bcr_status[test.idx]
  
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
  
  tr <- paste0('[Train] ', training.risk.score)
  names(tr) <- names(training.risk.score)
  
  test <- paste0('[Test] ', test.risk.score)
  names(test) <- names(test.risk.score)
  
  cv[[k]] <- c(tr, test)[samples]
  
  riskScore <- c(riskScore, test.risk.score)
  
  trainingRiskGroup <- c(trainingRiskGroup, training.risk.group)
  testRiskGroup <- c(testRiskGroup, test.risk.group)
  
}

cv <- do.call(cbind, cv)
colnames(cv) <- paste0('CV',1:10)
cv

riskScore <- riskScore[samples]
trainingRiskGroup <- trainingRiskGroup[samples]
testRiskGroup <- testRiskGroup[samples]

riskThreshold <- median(riskScore, na.rm=T)
riskGroup <- riskScore > riskThreshold

res <- data.frame(risk.score=riskScore, risk.group=riskGroup, training.risk.group=trainingRiskGroup, 
                  test.risk.group=testRiskGroup, bcr.time=pheno$time_to_bcr, bcr.status=pheno$bcr_status, cv)


##################################### survivalsvm

riskScore <- c()
trainingRiskGroup <- c()
testRiskGroup <- c()

cv <- list()

for (k in 1:nfold) {
  
  training.idx <- which(foldid!=k)
  test.idx <- which(foldid==k)
  
  training.geno <- geno[training.idx,genes]
  test.geno <- geno[test.idx,genes]
  
  filter <- which(apply(training.geno, 2, function(v) sum(v==mean(v))==length(v)))
  #filter <- which(apply(training.geno, 2, sd)==0)
  
  if (length(filter)>0) {
    training.geno <- training.geno[,-filter]
    test.geno <- test.geno[,-filter]
  }
  
  training.bcr.time <- pheno$time_to_bcr[training.idx]
  training.bcr.status <-  pheno$bcr_status[training.idx]
  
  test.bcr.time <- pheno$time_to_bcr[test.idx]
  test.bcr.status <-  pheno$bcr_status[test.idx]
  
  training.surv.data <- data.frame(training.geno, bcr.time=training.bcr.time, bcr.status=training.bcr.status)
  
  set.seed(777)
  svmfit <- survivalsvm(Surv(bcr.time, bcr.status) ~ ., data = training.surv.data, 
                        gamma.mu = 0.1, kernel = 'lin_kernel')
  
  training.risk.score <- predict(svmfit, data.frame(training.geno))
  training.risk.score
  training.risk.score <- training.risk.score$predicted[1,]
  
  training.risk.threshold <- median(training.risk.score, na.rm = T)
  training.risk.threshold
  
  test.risk.score <- predict(svmfit, data.frame(test.geno)) ###
  test.risk.score
  test.risk.score <- test.risk.score$predicted[1,]
  
  test.risk.threshold <- median(test.risk.score, na.rm = T)
  test.risk.threshold
  
  
  training.risk.group <- test.risk.score > training.risk.threshold
  training.risk.group
  
  test.risk.group <- test.risk.score > test.risk.threshold
  test.risk.group
  
  tr <- paste0('[Train] ', training.risk.score)
  names(tr) <- names(training.risk.score)
  
  test <- paste0('[Test] ', test.risk.score)
  names(test) <- names(test.risk.score)
  
  cv[[k]] <- c(tr, test)[samples]
  
  riskScore <- c(riskScore, test.risk.score)
  
  trainingRiskGroup <- c(trainingRiskGroup, training.risk.group)
  testRiskGroup <- c(testRiskGroup, test.risk.group)
  
}

cv <- do.call(cbind, cv)
colnames(cv) <- paste0('CV',1:10)
cv

riskScore <- riskScore[samples]
trainingRiskGroup <- trainingRiskGroup[samples]
testRiskGroup <- testRiskGroup[samples]

riskThreshold <- median(riskScore, na.rm=T)
riskGroup <- riskScore > riskThreshold

res <- data.frame(risk.score=riskScore, risk.group=riskGroup, training.risk.group=trainingRiskGroup, 
                  test.risk.group=testRiskGroup, bcr.time=pheno$time_to_bcr, bcr.status=pheno$bcr_status, cv)


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
