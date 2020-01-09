
setwd('~/bigdata/PCa/')
source('script/PCa_Functions.R')

library(Biobase)
library(caret)
library(readxl)
library(pROC)
library(randomForest)
library(plyr)
library(xgboost)
library(glmnet)
library(MASS) # lda
library(pls)
library(kernlab)

####### GSE46691#######

dataset <- 'GSE46691'

eSet <- readRDS(paste0('data/Database/Primary/', dataset, '_eSet.RDS'))
exprData <- exprs(eSet)
phenoData <- pData(eSet)
#View(phenoData)

#phenoData$bcr_status_3yr <- getEventFun(n=3, time.to.event=phenoData$time_to_bcr, status=phenoData$bcr_status)
#phenoData$bcr_status_5yr <- getEventFun(n=5, time.to.event=phenoData$time_to_bcr, status=phenoData$bcr_status)
#phenoData$bcr_status_10yr <- getEventFun(n=10, time.to.event=phenoData$time_to_bcr, status=phenoData$bcr_status)

#table(phenoData$bcr_status_3yr)
#table(phenoData$bcr_status_5yr)
#table(phenoData$bcr_status_10yr)
table(phenoData$metastasis_status)

#filter <- which(is.na(phenoData$bcr_status))
#filter
#phenoData <- phenoData[-filter,]
#exprData <- exprData[,-filter]

pheno <- as.matrix(paste0('E',phenoData$metastasis_status), drop=FALSE) ### 5 Year
pheno
geno <- as.matrix(t(exprData))
#geno <- scale(geno)

set.seed(777)
foldid <- generateCVFold(nsam=length(pheno), nfold=10)
#saveRDS(foldid, file=paste0('data/foldid/Foldid_', dataset, '_bcr_5yr.RDS'))


### Test a single signature
signature <- read_xlsx(path = 'data/Classifiers.xlsx', sheet='Oncotype')
signature <- signature$Ensembl

pred.test <- runCaretCV(geno = geno, pheno = pheno, foldid = foldid, signature = signature, model = 'glmnet')

table(pred=pred.test$ypred, true=pred.test$yobs)

roc.test <- roc(pred.test$yobs, pred.test$yprob, plot=FALSE, ci=TRUE, auc=TRUE)
auc <- roc.test$ci[2]
auc.ci.lower95 <- roc.test$ci[1]
auc.ci.upper95 <- roc.test$ci[3]
auc
auc.ci.lower95
auc.ci.upper95


# Published signatures
signatures <- c('Agell','Bibikova','Bismar','Decipher','Ding','Glinsky','Irshad',
                'Jennifer','Jia','Kamoun','Long','Luca','Mo','Nakagawa','Olmos',
                'Oncotype','Penney','Planche','Prolaris','Ramaswamy','Ramos_Montoya',
                'Ross_Adams','Ross_Robert','Sharma','Talantov','Varambally','Wu','Yang',
                'Yu')

# Random signatures
random.signatures <- readRDS(file='data/Random_Signatures.RDS')

models <- c('glmnet','svmLinear','svmRadial','svmPoly','rf','pls','lda','xgbLinear', 'xgbTree') # dnn
models <- models[-1]
models

for (model in models) {
  message(model)
  
  ### published signatures
  for (signature.name in signatures) {
    message(signature.name)
    
    signature <- read_xlsx(path = 'data/Classifiers.xlsx', sheet=signature.name)
    signature.genes <- signature$Ensembl
    
    pred.test <- runCaretCV(geno = geno, pheno = pheno, foldid = foldid, signature = signature.genes, model = model)
    write.csv(pred.test, file=paste0('report/Classification/CV10/YR5/Classification_', model, '_', signature.name, '_', dataset, '.csv'),
              quote = F,)
  }
  
  ### random signatures
  for (signature.name in names(random.signatures)) {
    message (signature.name)
    
    signature.genes <- random.signatures[[signature.name]]
    
    pred.test <- runCaretCV(geno = geno, pheno = pheno, foldid = foldid, signature = signature.genes, model = model)
    write.csv(pred.test, file=paste0('report/Classification/CV10/YR5/Classification_', model, '_', signature.name, '_', dataset, '.csv'),
              quote = F)
  }
  
}


