######################################

gtf <- readRDS(file='data/GENCODE_Annotation_Human_V32.RDS')

classifier <- read_xlsx(path = 'data/Classifiers.xlsx', sheet='Oncotype')
classifier

classifier <- gtf[match(classifier$TCGA, gtf$gene_name),]
classifier

which(is.na(classifier$ensembl))

write.csv(classifier, file='data/Classifiers/Ross_Adams_GENCODE32.csv', quote=F)

gtf[gtf$ensembl=='ENSG00000182253',]

lucaClassifiers <- read_xlsx(path = 'data/Classifiers/Luca_2017_CancerMap_Suppl_Table2_Classifiers.xlsx')
unique(lucaClassifiers$Signature)
