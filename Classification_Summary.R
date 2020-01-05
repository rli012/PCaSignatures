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
      
      a <- sum(pred.test$yobs=='E1' & pred.test$ypred=='E1')
      b <- sum(pred.test$yobs=='E0' & pred.test$ypred=='E1')
      c <- sum(pred.test$yobs=='E1' & pred.test$ypred=='E0')
      d <- sum(pred.test$yobs=='E0' & pred.test$ypred=='E0')
      
      sensitivity <- a/(a+c)
      specificity <- d/(b+d)
      
      accuracy <- (a+d)/(a+b+c+d)
      
      roc.test <- roc(pred.test$yobs, pred.test$yprob, plot=FALSE, ci=TRUE, auc=TRUE)
      auc <- round(roc.test$ci[2],3)
      auc.ci.lower95 <- round(roc.test$ci[1],3)
      auc.ci.upper95 <- round(roc.test$ci[3],3)
      
      print (paste0(auc, ' (', auc.ci.lower95, '-', auc.ci.upper95, ')'))
      
      sumpred <- rbind(sumpred, c(model, dataset, signature.name, 
                                  accuracy, sensitivity, specificity,
                                  auc, auc.ci.lower95, auc.ci.upper95))
      
    }
    
    for (signature.name in names(random.signatures)) {
      
      message(signature.name)
      fl <- paste0('report/Classification/CV10/YR5/Classification_', model, '_', signature.name, '_', dataset, '.csv')
      
      if (!file.exists(fl)) {
        next
      }
      
      pred.test <- read.csv(fl, row.names=1)
      #print(table(pred=pred.test$ypred, true=pred.test$yobs))
      
      a <- sum(pred.test$yobs=='E1' & pred.test$ypred=='E1')
      b <- sum(pred.test$yobs=='E0' & pred.test$ypred=='E1')
      c <- sum(pred.test$yobs=='E1' & pred.test$ypred=='E0')
      d <- sum(pred.test$yobs=='E0' & pred.test$ypred=='E0')
      
      sensitivity <- a/(a+c)
      specificity <- d/(b+d)
      
      accuracy <- (a+d)/(a+b+c+d)
      
      roc.test <- roc(pred.test$yobs, pred.test$yprob, plot=FALSE, ci=TRUE, auc=TRUE)
      auc <- round(roc.test$ci[2],3)
      auc.ci.lower95 <- round(roc.test$ci[1],3)
      auc.ci.upper95 <- round(roc.test$ci[3],3)
      
      print (paste0(auc, ' (', auc.ci.lower95, '-', auc.ci.upper95, ')'))
      
      sumpred <- rbind(sumpred, c(model, dataset, signature.name, 
                                  accuracy, sensitivity, specificity,
                                  auc, auc.ci.lower95, auc.ci.upper95))
      
    }
    
  }
  
}

sumpred <- data.frame(sumpred, stringsAsFactors = F)
colnames(sumpred) <- c('model','dataset','signature',
                       'accuracy','sensitivity','specificity',
                       'auc','auc.lower95','auc.upper95')

sumpred[,4:9] <- apply(sumpred[,4:9], 2, as.numeric)
saveRDS(sumpred, 'report/sumpred_for_test.RDS')
