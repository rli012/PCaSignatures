
####################################################################
###      Evaluation of Survial Analysis Models (10-fold CV)      ###
####################################################################

setwd('~/bigdata/PCa/')
source('script/PCa_Functions.R')
source('~/bigdata/LABDATA/BLUPHAT/script/Helper_Functions.R')


# BRB-ArrayTools: https://brb.nci.nih.gov/BRB-ArrayTools/download.html
# randomForestSRC: https://kogalur.github.io/randomForestSRC/theory.html
# superpc: http://statweb.stanford.edu/~tibs/superpc/tutorial.html
# plsRcox: https://fbertran.github.io/plsRcox/
# BhGLM: https://github.com/nyiuab/BhGLM

#library(remotes)
#install_github("nyiuab/BhGLM", force=T, build_vignettes=T)
#BiocManager::install('survcomp')

library(Biobase)
library(readxl)
library(survminer)
library(survcomp)

library(survival) # coxph
library(glmnet) # coxnet
library(superpc) # pca
#library(survivalsvm) # svm
#library(randomForestSRC)
#library(plsRcox)

#library(BhGLM)

# Elastic Net/Lasso-Cox
# survivalSVM
# randomForestSRC
# plsRcox
# Superpc

# Bayesian Losso
# Neural Network Cox
# Dr. Yi


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


datasets <- c('TCGA_PRAD','GSE107299','GSE21034','DKFZ2018','GSE54460','GSE70768','GSE70769',
              'GSE94767','E_MTAB_6128','GSE116918_BCR') # , 'GSE116918_Metastasis'

dataset <- datasets[10]  #####

eSet <- readRDS(paste0('data/Database/Primary/', dataset, '_eSet.RDS'))
exprData <- exprs(eSet)
phenoData <- pData(eSet)
#View(phenoData)
table(phenoData$sample_type)

keep <- which(phenoData$sample_type=='Primary' | phenoData$sample_type=='Tumor' | phenoData$sample_type=='Tumour' )
exprData <- exprData[,keep]
phenoData <- phenoData[keep,]


filter <- which(duplicated(phenoData$patient_id, incomparables = NA))

if (length(filter)>0) {
  exprData <- exprData[,-filter]
  phenoData <- phenoData[-filter,]
}


#pheno$time_to_bcr <- pheno$time_to_metastasis
#pheno$bcr_status <- pheno$metastasis_status

filter <- which(is.na(phenoData$bcr_status) | is.na(phenoData$time_to_bcr))

if (length(filter)>0) {
  exprData <- exprData[,-filter]
  phenoData <- phenoData[-filter,]
}


pheno <- phenoData
dim(pheno)
pheno$time_to_bcr[pheno$time_to_bcr<=0] <- 0.01

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
                'Jennifer','Jia','Kamoun','Long','Luca','Mo','Nakagawa','Olmos',
                'Oncotype','Penney','Planche','Prolaris','Ramaswamy','Ramos_Montoya',
                'Ross_Adams','Ross_Robert','Sharma','Talantov','Varambally','Wu','Yang',
                'Yu')

# signatures <- c('Penney', 'Irshad', 'Ramos_Montoya','Ross_Adams','Sharma', 'Yang')  # n > p in most of the datasets

genes
for (signature.name in signatures) {
  message(signature.name)
  
  signature <- read_xlsx(path = 'data/Classifiers.xlsx', sheet=signature.name)
  signature.genes <- signature$Ensembl

  genes <- intersect(signature.genes, colnames(geno))

  ################################ CoxPH
  
  model <- 'CoxPH'
  print (model)

  riskScore <- c()
  trainingRiskGroup <- c()
  testRiskGroup <- c()

  cv <- list()

  for (k in 1:nfold) {

    training.idx <- which(foldid!=k)
    test.idx <- which(foldid==k)

    training.geno <- geno[training.idx,genes]
    test.geno <- geno[test.idx,genes]

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

    training.risk.score <- predict(coxtest, data.frame(training.geno), type="lp", se.fit=FALSE) # type='risk'
    training.risk.score

    training.risk.threshold <- median(training.risk.score, na.rm = T)
    training.risk.threshold

    test.risk.score <- predict(coxtest, data.frame(test.geno), type="lp", se.fit=FALSE)
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


  write.csv(res, file=paste0('report/Survival/CV10Scale/Signature/Survival_', model, '_', signature.name, '_', dataset, '.csv'),
            quote = F)

  
  ############################## Elastic Net (CoxNet)
  
  model <- 'CoxNet'
  print (model)
  
  riskScore <- c()
  trainingRiskGroup <- c()
  testRiskGroup <- c()
  
  cv <- list()
  
  for (k in 1:nfold) {
    
    training.idx <- which(foldid!=k)
    test.idx <- which(foldid==k)
    
    training.geno <- geno[training.idx,genes]
    test.geno <- geno[test.idx,genes]
    
    training.bcr.time <- pheno$time_to_bcr[training.idx]
    training.bcr.status <-  pheno$bcr_status[training.idx]
    
    test.bcr.time <- pheno$time_to_bcr[test.idx]
    test.bcr.status <-  pheno$bcr_status[test.idx]
    
    training.surv.data <- data.frame(training.geno, bcr.time=training.bcr.time, bcr.status=training.bcr.status)
    
    set.seed(777)
    cv.fit <- cv.glmnet(training.geno, Surv(training.bcr.time, training.bcr.status), family="cox", maxit = 1000)
    
    coeffs <- coef(cv.fit, s = cv.fit$lambda.min)
    coeffs <- as.numeric(coeffs)
    
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
  
  write.csv(res, file=paste0('report/Survival/CV10Scale/Signature/Survival_', model, '_', signature.name, '_', dataset, '.csv'),
            quote = F)
  
  ##################################### superpc
  
  model <- 'SuperPC'
  print (model)
  
  riskScore <- c()
  trainingRiskGroup <- c()
  testRiskGroup <- c()
  
  cv <- list()
  
  for (k in 1:nfold) {
    
    training.idx <- which(foldid!=k)
    test.idx <- which(foldid==k)
    
    training.geno <- geno[training.idx,genes]
    test.geno <- geno[test.idx,genes]
    
    training.bcr.time <- pheno$time_to_bcr[training.idx]
    training.bcr.status <-  pheno$bcr_status[training.idx]
    
    test.bcr.time <- pheno$time_to_bcr[test.idx]
    test.bcr.status <-  pheno$bcr_status[test.idx]
    
    training.data.pca <-list(x=t(training.geno), y=training.bcr.time, censoring.status= training.bcr.status, featurenames=genes)
    test.data.pca <- list(x=t(test.geno), y=test.bcr.time, censoring.status= test.bcr.status, featurenames=genes)
    
    set.seed(777)
    pca.fit <- superpc.train(training.data.pca, type="survival")
    
    training.risk.score <- superpc.predict(pca.fit, training.data.pca, training.data.pca, 
                                           threshold=0.3, n.components=1, prediction.type="continuous")
    
    training.risk.score <- training.risk.score$v.pred[,1]
    names(training.risk.score) <- rownames(training.geno)
    
    training.risk.threshold <- median(training.risk.score, na.rm = T)
    
    test.risk.score <- superpc.predict(pca.fit, training.data.pca, test.data.pca, 
                                       threshold=0.3, n.components=1, prediction.type="continuous")
    
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
  
  write.csv(res, file=paste0('report/Survival/CV10Scale/Signature/Survival_', model, '_', signature.name, '_', dataset, '.csv'),
            quote = F)
  
  ##################################### survivalsvm
  
  # Ctrl + Shift + C, comment out a block of code
  
  # model <- 'SurvivalSVM'
  # 
  # riskScore <- c()
  # trainingRiskGroup <- c()
  # testRiskGroup <- c()
  # 
  # cv <- list()
  # 
  # for (k in 1:nfold) {
  #   
  #   training.idx <- which(foldid!=k)
  #   test.idx <- which(foldid==k)
  #   
  #   training.geno <- geno[training.idx,genes]
  #   test.geno <- geno[test.idx,genes]
  #   
  #   training.bcr.time <- pheno$time_to_bcr[training.idx]
  #   training.bcr.status <-  pheno$bcr_status[training.idx]
  #   
  #   test.bcr.time <- pheno$time_to_bcr[test.idx]
  #   test.bcr.status <-  pheno$bcr_status[test.idx]
  #   
  #   training.surv.data <- data.frame(training.geno, bcr.time=training.bcr.time, bcr.status=training.bcr.status)
  #   
  #   set.seed(777)
  #   svmfit <- survivalsvm(Surv(bcr.time, bcr.status) ~ ., data = training.surv.data, 
  #                         gamma.mu = 0.1, kernel = 'lin_kernel')
  #   
  #   training.risk.score <- predict(svmfit, data.frame(training.geno))
  #   training.risk.score
  #   training.risk.score <- training.risk.score$predicted[1,]
  #   
  #   training.risk.threshold <- median(training.risk.score, na.rm = T)
  #   training.risk.threshold
  #   
  #   test.risk.score <- predict(svmfit, data.frame(test.geno)) ###
  #   test.risk.score
  #   test.risk.score <- test.risk.score$predicted[1,]
  #   
  #   test.risk.threshold <- median(test.risk.score, na.rm = T)
  #   test.risk.threshold
  #   
  #   
  #   training.risk.group <- test.risk.score > training.risk.threshold
  #   training.risk.group
  #   
  #   test.risk.group <- test.risk.score > test.risk.threshold
  #   test.risk.group
  #   
  #   tr <- paste0('[Train] ', training.risk.score)
  #   names(tr) <- names(training.risk.score)
  #   
  #   test <- paste0('[Test] ', test.risk.score)
  #   names(test) <- names(test.risk.score)
  #   
  #   cv[[k]] <- c(tr, test)[samples]
  #   
  #   riskScore <- c(riskScore, test.risk.score)
  #   
  #   trainingRiskGroup <- c(trainingRiskGroup, training.risk.group)
  #   testRiskGroup <- c(testRiskGroup, test.risk.group)
  #   
  # }
  # 
  # cv <- do.call(cbind, cv)
  # colnames(cv) <- paste0('CV',1:10)
  # cv
  # 
  # riskScore <- riskScore[samples]
  # trainingRiskGroup <- trainingRiskGroup[samples]
  # testRiskGroup <- testRiskGroup[samples]
  # 
  # riskThreshold <- median(riskScore, na.rm=T)
  # riskGroup <- riskScore > riskThreshold
  # 
  # res <- data.frame(risk.score=riskScore, risk.group=riskGroup, training.risk.group=trainingRiskGroup, 
  #                   test.risk.group=testRiskGroup, bcr.time=pheno$time_to_bcr, bcr.status=pheno$bcr_status, cv)
  # 
  # 
  # write.csv(res, file=paste0('report/Survival/CV10Scale/Signature/Survival_', model, '_', signature.name, '_', dataset, '.csv'),
  #           quote = F)
  
  
}
