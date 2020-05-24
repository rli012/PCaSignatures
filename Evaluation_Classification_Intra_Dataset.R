
#######################################################################
###       Evaluation of Classification Models (Intra-Dataset)       ###
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
              'GSE94767','E_MTAB_6128','GSE116918_BCR') # , 'GSE116918_Metastasis'

dataset <- datasets[1]  #####

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

# Random signatures
#random.signatures <- readRDS(file='data/Random_Signatures.RDS')

models <- c('glmnet','svmLinear','svmRadial','svmPoly','rf','pls','lda','xgbLinear','xgbTree') # dnn

for (model in models) {
  message(model)
  
  ### published signatures
  for (signature.name in signatures) {
    message(signature.name)
    
    signature <- read_xlsx(path = 'data/PCa_Prognosis_Signatures.xlsx', sheet=signature.name)
    signature.genes <- signature$Ensembl
    
    pred.test <- runCaretCV(geno = geno, pheno = pheno, foldid = foldid, signature = signature.genes, model = model)
    write.csv(pred.test, file=paste0('report/Classification/CV10Scale/Signature/Classification_', model, '_', signature.name, '_', dataset, '.csv'),
              quote = F,)
  }
  
  # ### random signatures
  # for (signature.name in names(random.signatures)) {
  #   message (signature.name)
  #   
  #   signature.genes <- random.signatures[[signature.name]]
  #   
  #   pred.test <- runCaretCV(geno = geno, pheno = pheno, foldid = foldid, signature = signature.genes, model = model)
  #   write.csv(pred.test, file=paste0('report/Classification/CV10Scale/Signature/Classification_', model, '_', signature.name, '_', dataset, '.csv'),
  #             quote = F)
  # }
  
}
