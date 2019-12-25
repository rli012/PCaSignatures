
getGENCODEAnnotation <- function(species='human', release='31', type='gene') {
  
  # species: human, mouse
  # release: human 31, mouse M20
  # type: gene, transcript
  
  baseurl <- 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/'
  
  gtf.file <- paste0(baseurl, 'Gencode_', species, '/release_', release, '/gencode.v', release, '.annotation.gtf.gz')
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


getSurvivalFun <- function(n=3, days.survival, days.followup) {
  
  idx <- which(days.followup > n*365)
  
  if (length(idx) > 0) {
    days.survival[idx] <- days.followup[idx]
  } 
  
  days.survival <- ifelse(days.survival > n*365, 'E0', 'E1')
  
  return(days.survival)
}

#clinicalDa$BCR_Status_5YR <- getSurvivalFun(n=5, 
#                                            days.survival = clinicalDa$days_to_first_biochemical_recurrence, 
#                                            days.followup = clinicalDa$days_to_last_followup)


generateCVFold <- function(nsam, nfold=10) {
  sample(rep(1:nfold,ceiling(nsam/nfold))[1:nsam])
}

#foldid <- generateCVFold(nsam=153, nfold=10)
