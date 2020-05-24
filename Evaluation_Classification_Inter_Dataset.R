``
#######################################################################
###       Evaluation of Classification Models (Inter-Dataset)       ###
#######################################################################

setwd('~/bigdata/PCa/')
source('script/Helper_Functions.R')

library(Biobase)
library(caret)
library(readxl)
library(randomForest)
library(xgboost)
library(glmnet)
library(MASS) # lda
library(pls)
library(kernlab)

options(stringsAsFactors = F)

#############################################################################

datasets <- c('TCGA_PRAD','GSE107299','GSE21034','DKFZ2018','GSE54460','GSE70768','GSE70769',
                  'GSE94767','E_MTAB_6128','GSE116918_BCR')

training.set <- datasets[1]
test.sets <- datasets[which(!datasets %in% training.set)]

eSet <- readRDS(paste0('data/Database/Primary/', training.set, '_eSet.RDS'))
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

# time.bcr.set <- c('TCGA_PRAD','GSE107299','GSE21034','DKFZ2018','GSE54460','GSE70768','GSE70769',
#                   'GSE94767','E_MTAB_6128','GSE116918_BCR')
# time.met.set <- c('GSE116918_Metastasis')
# 
# binary.bcr.set <- c('GSE25136','GSE41408_BCR')
# binary.met.set <- c('GSE41408_Metastasis','GSE46691','GSE51066')
# 
# special.set <- 'GSE37199'
# 
# test.sets <- c(time.bcr.set, time.met.set, binary.bcr.set, binary.met.set, special.set)
# test.sets <- test.sets[which(!test.sets %in% training.set)]
# test.sets

# if (training.set %in% time.bcr.set) {
#   #phenoData$bcr_status_3yr <- getEventFun(n=3, time.to.event=phenoData$time_to_bcr, status=phenoData$bcr_status)
#   phenoData$risk_group <- getEventFun(n=5, time.to.event=phenoData$time_to_bcr, status=phenoData$bcr_status)
#   #phenoData$bcr_status_10yr <- getEventFun(n=10, time.to.event=phenoData$time_to_bcr, status=phenoData$bcr_status)
#   
# } else if (training.set %in% time.met.set) {
#   phenoData$risk_group <- getEventFun(n=5, time.to.event=phenoData$time_to_metastasis, status=phenoData$metastasis_status)
#   
# } else if (training.set %in% binary.bcr.set) {
#   phenoData$risk_group <- ifelse(is.na(phenoData$bcr_status), NA, paste0('E',phenoData$bcr_status))
#   
# } else if (training.set %in% binary.met.set) {
#   phenoData$risk_group <- ifelse(is.na(phenoData$metastasis_status), NA, paste0('E',phenoData$metastasis_status))
#   
# } else if (training.set %in% special.set) {
#   phenoData$risk_group <- ifelse(phenoData$risk_group=='Good prognosis', 'E0', 'E1')
# }

phenoData$risk_group <- getEventFun(n=5, time.to.event=phenoData$time_to_bcr, status=phenoData$bcr_status)

keep <- which(!is.na(phenoData$risk_group))
phenoData <- phenoData[keep,]
exprData <- exprData[,keep]
#phenoData

training.pheno <- as.matrix(phenoData$risk_group, drop=FALSE) ### 5 Year
training.geno.all <- as.matrix(t(exprData))
training.geno.all <- scale(training.geno.all)

#set.seed(777)
#foldid <- generateCVFold(n=length(pheno), nfold=10)
#saveRDS(foldid, file=paste0('data/foldid/Foldid_', dataset, '_bcr_5yr.RDS'))


#############################################################################

# Published signatures
signatures <- c('Agell','Bibikova','Bismar','Decipher','Ding','Glinsky','Irshad',
                'Jennifer','Jia','Kamoun','Li','Long','Luca','Mo','Nakagawa','Olmos',
                'Oncotype','Penney','Planche','Prolaris','Ramaswamy','Ramos_Montoya',
                'Ross_Adams','Ross_Robert','Sharma','Talantov','Varambally','Wu','Yang',
                'Yu')

# Random signatures
# random.signatures <- readRDS(file='data/Random_Signatures.RDS')

models <- c('glmnet','svmLinear','svmRadial','svmPoly','rf','pls','lda','xgbLinear', 'xgbTree') # dnn

for (model in models) {
  message(model)
  
  ### published signatures
  for (signature.name in signatures) {
    message(signature.name)
    
    signature <- read_xlsx(path = 'data/PCa_Prognosis_Signatures.xlsx', sheet=signature.name)
    signature.genes <- signature$Ensembl
    
    for (test.set in test.sets) {
      eSet <- readRDS(paste0('data/Database/Primary/', test.set, '_eSet.RDS'))
      
      exprData <- exprs(eSet)
      phenoData <- pData(eSet)
      keep <- which(phenoData$sample_type=='Primary' | phenoData$sample_type=='Tumor' | phenoData$sample_type=='Tumour' )
      
      exprData <- exprData[,keep]
      phenoData <- phenoData[keep,]
      
      filter <- which(duplicated(phenoData$patient_id, incomparables = NA))
      
      if (length(filter)>0) {
        exprData <- exprData[,-filter]
        phenoData <- phenoData[-filter,]
      }
      
      # if (test.set %in% time.bcr.set) {
      #   #phenoData$bcr_status_3yr <- getEventFun(n=3, time.to.event=phenoData$time_to_bcr, status=phenoData$bcr_status)
      #   phenoData$risk_group <- getEventFun(n=5, time.to.event=phenoData$time_to_bcr, status=phenoData$bcr_status)
      #   #phenoData$bcr_status_10yr <- getEventFun(n=10, time.to.event=phenoData$time_to_bcr, status=phenoData$bcr_status)
      #   
      # } else if (test.set %in% time.met.set) {
      #   phenoData$risk_group <- getEventFun(n=5, time.to.event=phenoData$time_to_metastasis, status=phenoData$metastasis_status)
      #   
      # } else if (test.set %in% binary.bcr.set) {
      #   phenoData$risk_group <- paste0('E',phenoData$bcr_status) ### 5 Year
      #   
      # } else if (test.set %in% binary.met.set) {
      #   phenoData$risk_group <- paste0('E',phenoData$metastasis_status) ### 5 Year
      #   
      # } else if (test.set %in% special.set) {
      #   phenoData$risk_group <- ifelse(phenoData$risk_group=='Good prognosis', 'E0', 'E1')
      # }
      
      phenoData$risk_group <- getEventFun(n=5, time.to.event=phenoData$time_to_bcr, status=phenoData$bcr_status)
      
      keep <- which(!is.na(phenoData$risk_group))
      phenoData <- phenoData[keep,]
      exprData <- exprData[,keep]
      
      test.pheno <- as.matrix(phenoData$risk_group, drop=FALSE) ### 5 Year
      test.geno <- as.matrix(t(exprData))
      test.geno <- scale(test.geno)
      
      genes <- Reduce(intersect, list(signature.genes, colnames(training.geno.all), colnames(test.geno)))
      
      training.geno <- training.geno.all[,genes]
      test.geno <- test.geno[,genes]
      
      filter <- which(apply(training.geno, 2, function(v) sum(v==mean(v))==length(v)))
      #filter <- which(apply(training.geno, 2, sd)==0)
      
      if (length(filter)>0) {
        training.geno <- training.geno[,-filter]
        test.geno <- test.geno[,-filter]
      }
      
      #training.geno <- scale(training.geno)
      #test.geno <- scale(test.geno)
      
      tr.control <- trainControl(
        method = 'cv', #"repeatedcv",
        number = 10,
        #repeats = 3, 
        classProbs = TRUE,
        #returnResamp="all",
        search = 'grid', # random
        summaryFunction=twoClassSummary)
      
      set.seed(777)
      
      if (model %in% c('glmnet','svmLinear','svmRadial','rf','pls','lda')) {
        model.fit <- train(training.geno, training.pheno,
                           #y1 ~ ., data = tr.data,
                           method = model,
                           metric = 'ROC',
                           #preProc = c("center", "scale"),
                           #tuneGrid = tune.grid,
                           tuneLength = 10, 
                           trControl = tr.control)
        
      } else if (model %in% c('svmPoly','xgbLinear', 'xgbTree','dnn')) {
        model.fit <- train(training.geno, training.pheno,
                           #y1 ~ ., data = tr.data,
                           method = model,
                           metric = 'ROC',
                           #preProc = c("center", "scale"),
                           #tuneGrid = tune.grid,
                           #tuneLength = 10, 
                           trControl = tr.control)
      }
      
      yprob <- predict(model.fit, newdata = test.geno, type = "prob")
      yprob <- yprob$E1
      ypred <- predict(model.fit, newdata = test.geno, type = "raw")
      
      pred.test <- data.frame(yobs=test.pheno, ypred=ypred, yprob, row.names=rownames(phenoData), stringsAsFactors = F)
      write.csv(pred.test, file=paste0('report/Classification/TrainTestScale/Signature/Classification_', 
                                       model, '_', signature.name, '_', training.set, '_', test.set, '_TrainTest_Signature.csv'),
                quote = F)
    }
  }
}
