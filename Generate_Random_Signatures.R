
### Generate random signatures

datasets <- list.files('data/Database/Primary/')

gene.list <- c()
gene.list.small <- c()

for (dataset in datasets) {
  eSet <- readRDS(file.path('data/Database/Primary/', dataset))
  expr <- exprs(eSet)
  
  if (nrow(expr)<10000) {
    gene.list.small[[dataset]] <- rownames(expr)
    next
  }
  
  gene.list[[dataset]] <- rownames(expr)
}

genes <- Reduce(intersect, gene.list)
genes

random.signatures <- list()
random.signatures.wu <- list()

for (n in c(10,30,50,100)) {
  
  for (i in 1:3) {
    signature.name <- paste('Rand',n,i,sep='_')
    random.signatures[[signature.name]] <- sample(genes, n)
    
    signature.name <- paste('Wu_Rand',n,i,sep='_')
    random.signatures.wu[[signature.name]] <- sample(gene.list.small[[1]], n)
  }
}

saveRDS(random.signatures, file='data/Random_Signatures.RDS')
saveRDS(random.signatures.wu, file='data/Random_Signatures_Wu.RDS')


###
signatures <- c('Agell','Bibikova','Bismar','Decipher','Ding','Glinsky','Irshad',
                'Jennifer','Jia','Kamoun','Long','Luca','Mo','Nakagawa','Olmos',
                'Oncotype','Penney','Planche','Prolaris','Ramaswamy','Ramos_Montoya',
                'Ross_Adams','Ross_Robert','Sharma','Talantov','Varambally','Wu','Yang',
                'Yu')

signature.genes <- c()
for (signature.name in signatures) {
    message (signature.name)
    
    signature <- read_xlsx(path = 'data/Classifiers.xlsx', sheet=signature.name)
    signature.genes <- unique(c(signature.genes, signature$Ensembl))
    
}

signature.genes

ensembl <- readRDS('data/Annotation/ENSEMBL_Annotation_Human_V98_20191230.RDS')
ensembl

for (signature in names(random.signatures)) {
  idx <- match(random.signatures[[signature]], ensembl$ensembl_gene_id)
  #print (ensembl$external_gene_name[idx])
  print (sum(ensembl$external_gene_name[idx] %in% signature.genes))
}

for (signature in names(random.signatures.wu)) {
  idx <- match(random.signatures.wu[[signature]], ensembl$ensembl_gene_id)
  #print (ensembl$external_gene_name[idx])
  print (sum(ensembl$external_gene_name[idx] %in% signature.genes))
}


