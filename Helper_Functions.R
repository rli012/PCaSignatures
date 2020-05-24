###

runCaretCV <- function(geno, pheno, foldid, signature=NULL, model=NULL, seed=777) {
  nfold <- length(unique(foldid))
  
  if (is.null(signature)) {
    signature <- colnames(geno)
  } else {
    signature <- intersect(signature, colnames(geno))
  }
  
  samples <- NULL
  OUTPUT <- NULL
  
  for (k in 1:nfold) {
    message(paste(rep(c('+',k,'+'), each=20), collapse = '')) 
    i1<-which(foldid!=k)
    i2<-which(foldid==k)
    
    x1<-geno[i1,signature,drop=F]
    y1<-pheno[i1,,drop=F]
    
    x2<-geno[i2,signature,drop=F]
    y2<-pheno[i2,,drop=F]
    
    samples <- c(samples, rownames(geno)[i2])
    
    tr.control <- trainControl(
      method = 'cv', #"repeatedcv",
      number = 10,
      #repeats = 3, 
      classProbs = TRUE,
      #returnResamp="all",
      search = 'grid', # random
      summaryFunction=twoClassSummary)
    
    set.seed(seed)
    if (model %in% c('glmnet','svmLinear','svmRadial','rf','pls','lda')) {
      model.fit <- train(x1, y1,
                         #y1 ~ ., data = tr.data,
                         method = model,
                         metric = 'ROC',
                         #preProc = c("center", "scale"),
                         #tuneGrid = tune.grid,
                         tuneLength = 10, 
                         trControl = tr.control)
      
    } else if (model %in% c('svmPoly','xgbLinear', 'xgbTree','dnn')) {
      model.fit <- train(x1, y1,
                         #y1 ~ ., data = tr.data,
                         method = model,
                         metric = 'ROC',
                         #preProc = c("center", "scale"),
                         #tuneGrid = tune.grid,
                         #tuneLength = 10, 
                         trControl = tr.control)
      
    }
    
    yprob <- predict(model.fit, newdata = x2, type = "prob")
    yprob=yprob$E1
    ypred <- predict(model.fit, newdata = x2, type = "raw")
    
    OUTPUT <- rbind(OUTPUT, data.frame(yobs=y2, ypred=ypred, yprob))
  }
  
  OUTPUT <- data.frame(OUTPUT, stringsAsFactors = F)
  colnames(OUTPUT) <- c('yobs', 'ypred', 'yprob')
  rownames(OUTPUT) <- samples
  
  return (OUTPUT)
}



getGENCODEAnnotation <- function(species='human', release='32', type='gene', gtf.file=NULL) {
  
  # species: human, mouse
  # release: human 31, mouse M20
  # type: gene, transcript
  
  if (is.null(gtf.file)) {
    baseurl <- 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/'
    
    gtf.file <- paste0(baseurl, 'Gencode_', species, '/release_', release, '/gencode.v', release, '.annotation.gtf.gz')
  }
  
  gtf <- readGFF(gtf.file, version=2L)
  
  if (type!='all') {
    gtf <- gtf[gtf$type==type,]
    ensembl <- sapply(gtf$gene_id, function(x) strsplit(x, '.', fixed=T)[[1]][1])
    gtf <- add_column(gtf, ensembl, .before = 'gene_id')
  }
  
  return(gtf)
}

#gtf <- getGENCODEAnnotation(species='human', release='32', type='gene')
#rownames(gtf) <- ifelse(duplicated(gtf$ensembl), gtf$gene_id, gtf$ensembl)
#gtf <- add_column(.data = gtf, .after = 'end', length=gtf$end-gtf$start+1)
#saveRDS(gtf, file='data/GENCODE_Annotation_Human_V32.RDS')

#library(rtracklayer)
#gtf <- getGENCODEAnnotation(species='human', release='32', type='gene', gtf.file = 'data/Homo_sapiens.GRCh37.62.gtf.gz')



library(biomaRt)

getENSEMBLAnnotation <- function(attributes=NULL) {
  
  if (is.null(attributes)) {
    attributes <- c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol', 
                    'external_gene_name', 'description', 'gene_biotype')
  }
  
  #listMarts()
  
  ensembl=useMart("ensembl")
  #ensembl
  datasets <- listDatasets(ensembl)
  #head(datasets)
  
  ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
  #attributes <- ensembl@attributes
  #View(attributes)
  
  #affyids=c("202763_at","209310_s_at","207500_at")
  annotation <- getBM(attributes=attributes,
                      #filters = 'affy_hg_u133_plus_2', 
                      #values = affyids, 
                      mart = ensembl)
  
  return(annotation)
}

#ensembl <- getENSEMBLAnnotation()
#gtf <- readRDS(file='data/GENCODE_Annotation_Human_V32.RDS')
#sum(gtf$ensembl %in% ensembl$ensembl_gene_id)
#length(gtf$ensembl)
#saveRDS(ensembl, file='data/ENSEMBL_Annotation_Human_V98.RDS')

#attributes <- c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name', 'illumina_humanht_12_v4')
#ILLUMINA.HumanHT.12.V4 <- getENSEMBLAnnotation(attributes = attributes)


################

#gene.history <- read.table('data/gene_history.gz', header=T, stringsAsFactors = F, sep='\t', quote='', comment.char='')
#gene.history <- gene.history[gene.history$X.tax_id=='9606',]
#dim(gene.history)

#write.table(gene.history, file='data/Homo_sapiens. gene_history_20191226.txt', quote=F, sep='\t')



getEventFun <- function(n=3, time.to.event, status) {

  idx <- which(time.to.event < n*12 & status != 1)
  
  if (length(idx) > 0) {
    time.to.event[idx] <- NA
  } 
  
  group <- ifelse(time.to.event > n*12, 'E0', 'E1')
  
  return(group)
}

#clinicalDa$BCR_Status_5YR <- getEventFun(n=3, time.to.event=phenoData$time_to_bcr, status=phenoData$bcr_status)



generateCVFold <- function(n, nfold=10) {
  sample(rep(1:nfold,ceiling(n/nfold))[1:n])
}

#foldid <- generateCVFold(nsam=153, nfold=10)


### 10 repeats
generateCV <- function(n, nfold=10, repeats=10) {
  foldidID<-lapply(1:repeats,function(i){
    sample(rep(1:nfold,ceiling(n/nfold))[1:n])
  })
  return (foldidID)
}



### select most informative probe, MAX IQR

selectProbeFun <- function(expr) {
  expr$IQR <- apply(expr[,-which(colnames(expr)=='ID')], 1, IQR)
  
  expr <- expr %>% group_by(ID) %>% 
    filter(row_number() == which.max(IQR)) %>%
    column_to_rownames('ID')
  
  expr <- expr[,-which(colnames(expr)=='IQR')]
  
  return(expr)
  
}

#t <- expr %>% group_by(ID) %>% 
#  filter(row_number() == which.max(IQR)) %>%
#  column_to_rownames('ID')

#probe.idx <- expr %>% group_by(ID) %>%
#  summarise(idx=which.max(IQR))

#g <- rownames(t)[1]

#probe.idx[probe.idx$ID==g,]
#expr[expr$ID==g,1:5]
#t[1,1:5]
