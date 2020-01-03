
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


getEventFun <- function(n=3, time.to.event, status) {

  idx <- which(time.to.event < n*12 & status != 1)
  
  if (length(idx) > 0) {
    time.to.event[idx] <- NA
  } 
  
  group <- ifelse(time.to.event > n*12, 'E0', 'E1')
  
  return(group)
}

#clinicalDa$BCR_Status_5YR <- getEventFun(n=3, time.to.event=phenoData$time_to_bcr, status=phenoData$bcr_status)


generateCVFold <- function(nsam, nfold=10) {
  sample(rep(1:nfold,ceiling(nsam/nfold))[1:nsam])
}

#foldid <- generateCVFold(nsam=153, nfold=10)
