sumpred <- c()

for (model in models) {
  for (dataset in datasets) {
    
    message(dataset)
    
    for (signature.name in signatures) {
      
      message(signature.name)
      fl <- paste0('report/Classification/CV10/YR5/Classification_', model, '_', signature.name, '_', dataset, '.csv')
      
      if (!file.exists(fl)) {
        next
      }
      pred.test <- read.csv(fl, row.names=1)
      #print(table(pred=pred.test$ypred, true=pred.test$yobs))
      roc.test <- roc(pred.test$yobs, pred.test$yprob, plot=FALSE, ci=TRUE, auc=TRUE)
      auc <- round(roc.test$ci[2],3)
      auc.ci.lower95 <- round(roc.test$ci[1],3)
      auc.ci.upper95 <- round(roc.test$ci[3],3)
      
      print (paste0(auc, ' (', auc.ci.lower95, '-', auc.ci.upper95, ')'))
      
      sumpred <- rbind(sumpred, c(model, dataset, signature.name, auc, auc.ci.lower95, auc.ci.upper95))
      
    }
    
    for (signature.name in names(random.signatures)) {
      
      message(signature.name)
      fl <- paste0('report/Classification/CV10/YR5/Classification_', model, '_', signature.name, '_', dataset, '.csv')
      
      if (!file.exists(fl)) {
        next
      }
      
      pred.test <- read.csv(fl, row.names=1)
      #print(table(pred=pred.test$ypred, true=pred.test$yobs))
      roc.test <- roc(pred.test$yobs, pred.test$yprob, plot=FALSE, ci=TRUE, auc=TRUE)
      auc <- round(roc.test$ci[2],3)
      auc.ci.lower95 <- round(roc.test$ci[1],3)
      auc.ci.upper95 <- round(roc.test$ci[3],3)
      
      print (paste0(auc, ' (', auc.ci.lower95, '-', auc.ci.upper95, ')'))
      
      sumpred <- rbind(sumpred, c(model, dataset, signature.name, auc, auc.ci.lower95, auc.ci.upper95))
      
    }
    
  }
  
}

sumpred <- data.frame(sumpred, stringsAsFactors = F)
colnames(sumpred) <- c('model','dataset','signature','auc','auc.lower95','auc.upper95')
sumpred$auc <- as.numeric(sumpred$auc)
saveRDS(sumpred, 'report/sumpred_for_test.RDS')

