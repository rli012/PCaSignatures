
setwd('~/bigdata/PCa/')

library(GDCRNATools)
library(edgeR)
library(limma)
library(GEOquery)
library(oligo)
library(tibble) # column_to_rownames()
library(readxl)
library(stringr) # str_extract()


### Primary Datasets

# TCGA-PRAD     == Cell 2015                == GDC         == RNAseq: CPM              == Time to BCR
# CPC-GENE      == Nature_2017              == GSE107299   == HuGene 2.0 ST & HTA 2.0  == Time to BCR
# MSKCC         ==_Cancer_Cell_2010         == GSE21034    == HuEx 1.0 ST              == Time to BCR
# DKFZ          == Cancer Cell 2018         == cBioPortal  == RNAseq: RPKM             == Time to BCR
# GSE54460      == Cancer Research 2014     == GSE54460    == RNAseq: CPM              == Time to BCR
# CamCap        == Nature Genetics 2016     == GSE70768    == Illumina HumanHT-12 V4.0 == Time to BCR
# Stockholm     == Nature Genetics 2016     == GSE70769    == Illumina HumanHT-12 V4.0 == Time to BCR
# DESNT         == Eur Urol Focus 2017      == GSE94767    == HuEx 1.0 ST              == Time to BCR
# E-MTAB-6128   == Annuals of Oncology 2018 == E-MTAB-6128 == HuGene 2.0 ST            == Time to BCR
# Belfast       == Journal of Clinical Oncology 2017  == GSE116918 == ADXPCv1a520642   == Time to BCR/Metastasis

# GSE25136      == Prostate 2019            == GSE25136    == HG-U133A                 == BCR Status
# GSE44353      == PNAS 2013 == GSE44353    == Prostate Cancer DASL Panel 1.5K         == BCR Status (1500 genes) xxxx
# Rotterdam     == Int J Cancer 2013        == GSE41408    == HuEx 1.0 ST              == BCR/Metastasis Status
# Mayo Clinic I == Plos One 2013 (Decipher) == GSE46691    == HuEx 1.0 ST              == Metastasis Status
# GSE51066      == Prostate Cancer Prostatic Dis 2014 == GSE51066 == HuEx 1.0 ST       == Metastasis Status
# GSE37199      == Lancet Oncol 2012        == GSE37199    == HG-U133 Plus 2.0         == Prognosis Status (Good; Advanced Castration Resistant)


##########################################################################################
# ======================================== TCGA ======================================== #
##########################################################################################


#=========================================================================================
project <- 'TCGA-PRAD'
rnadir <- paste('data', project, 'RNAseq', sep='/')
mirdir <- paste('data', project, 'miRNAs', sep='/')

####### Download RNAseq data #######
gdcRNADownload(project.id     = project, 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = rnadir)


####### Download mature miRNA data #######
gdcRNADownload(project.id     = project, 
               data.type      = 'miRNAs', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = mirdir)



####### Download clinical data #######
clinicaldir <- paste('data', project, 'Clinical', sep='/')
gdcClinicalDownload(project.id     = project, 
                    write.manifest = FALSE,
                    method         = 'gdc-client',
                    directory      = clinicaldir)



####### Parse RNAseq metadata #######
metaMatrix.RNA <- gdcParseMetadata(project.id = project,
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)

####### Filter duplicated samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)

####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)

#saveRDS(metaMatrix.RNA, file='data/rData/Metadata_RNAseq_TCGA_PRAD.RDS')

####### Parse miRNAs metadata #######
metaMatrix.MIR <- gdcParseMetadata(project.id = project,
                                   data.type  = 'miRNAs', 
                                   write.meta = FALSE)

####### Filter duplicated samples in miRNAs metadata #######
metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)

####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in miRNAs metadata #######
metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)


####### Merge RNAseq data #######
rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA, 
                         path      = rnadir, # the folder in which the data stored
                         organized = FALSE, # if the data are in separate folders
                         data.type = 'RNAseq')

saveRDS(rnaCounts, file='data/TCGA-PRAD/RNAseq_Counts_TCGA_PRAD.RDS')

####### Merge miRNAs data #######
mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
                         path      = mirdir, # the folder in which the data stored
                         organized = FALSE, # if the data are in separate folders
                         data.type = 'miRNAs')

saveRDS(mirCounts, file='data/TCGA-PRAD/miRNA_Counts_TCGA_PRAD.RDS')



####### Merge clinical data #######
clinicalDa <- gdcClinicalMerge(path = clinicaldir, key.info = TRUE)
clinicalDa[1:6,5:10]

View(clinicalDa)


saveRDS(clinicalDa, file='data/TCGA-PRAD/Clinical_TCGA_PRAD_11262019.RDS')


sum(metaMatrix.MIR$sample_type=='PrimaryTumor')


pheno <- read.table('data/TCGA-PRAD/phenotype.TCGA-PRAD.final.txt', header=T, stringsAsFactors = F,
                    sep='\t', row.names = 1)

saveRDS(pheno, file='data/TCGA-PRAD/Clinical_TCGA_PRAD_With_PreopPSA_and_BCR.RDS')


##################

rnaCounts <- readRDS('data/TCGA-PRAD/RNAseq_Counts_TCGA_PRAD.RDS')
metaMatrix.RNA <- readRDS('data/TCGA-PRAD/Metadata_RNAseq_TCGA_PRAD.RDS')
mirCounts <- readRDS('data/TCGA-PRAD/miRNA_Counts_TCGA_PRAD.RDS')
metaMatrix.MIR <- readRDS('data/TCGA-PRAD/Metadata_miRNAs_TCGA_PRAD.RDS')
clinicalDa <- readRDS('data/TCGA-PRAD/Clinical_TCGA_PRAD_With_PreopPSA_and_BCR.RDS')

dge <-  DGEList(counts = rnaCounts)

### TMM normalization
dge = calcNormFactors(dge, method = 'TMM')

exprLogCPM <- edgeR::cpm(dge,log = TRUE) ### for visualization
dim(exprLogCPM)

saveRDS(exprLogCPM, file='data/rData/Expression_LogCPM_All_Genes_TCGA_PRAD.RDS')


### Filter out low-expression genes (cpm>1 in at least 50% of the samples)
keep <- rowSums(edgeR::cpm(dge) > 1) >= 0.5*ncol(rnaCounts)
sum(keep)
dge <- dge[keep,,keep.lib.sizes = TRUE]

### Voom normalization
#v <- voom(dge, design=NULL, plot = FALSE)

#exprAfterVoom <- v$E ### for visualization
exprLogCPMFilterLow <- edgeR::cpm(dge,log = TRUE) ### for visualization
exprLogCPMFilterLow

saveRDS(exprLogCPMFilterLow, file='data/rData/Expression_LogCPM_Filter_Low_Genes_TCGA_PRAD.RDS')



##########################################################################################
# ================================== Metadata from GEO ================================= #
##########################################################################################
###
gse <- 'GSE107299'

# GSE107299-GPL16686
seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE
seriesMatrix
# GSE107299-GPL16686
phenoData <- pData(seriesMatrix[[1]])
# GSE107299-GPL17586
phenoData <- pData(seriesMatrix[[2]])

View(phenoData)

phenoData <- phenoData[,c(1,2,8,36:38)]
colnames(phenoData)

colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
colnames(phenoData) <- gsub(' ', '_', colnames(phenoData))

colnames(phenoData)[3] <- 'tissue'

saveRDS(phenoData, file=paste0('data/rData/', gse, '_GPL16686_Sample_Information.RDS'))
saveRDS(phenoData, file=paste0('data/rData/', gse, '_GPL17586_Sample_Information.RDS'))


###
gse <- 'GSE21034'
seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE
names(seriesMatrix)
phenoData <- pData(seriesMatrix[[1]])
View(phenoData)
colnames(phenoData)

phenoData <- phenoData[,c(1,2,8,24,37:45)]


colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
colnames(phenoData) <- gsub(' ', '_', colnames(phenoData))

colnames(phenoData)[3] <- 'tissue'

colnames(phenoData)[8] <- 'clinical_t_stage'
colnames(phenoData)[10] <- 'pathological_t_stage'


saveRDS(phenoData, file=paste0('data/rData/', gse, '_GPL10264_Sample_Information.RDS'))


phenoData <- pData(seriesMatrix[[2]])
View(phenoData)
colnames(phenoData)

phenoData <- phenoData[,c(1,2,8,24,37:45)]


colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
colnames(phenoData) <- gsub(' ', '_', colnames(phenoData))

colnames(phenoData)[3] <- 'tissue'

colnames(phenoData)[8] <- 'clinical_t_stage'
colnames(phenoData)[10] <- 'pathological_t_stage'


saveRDS(phenoData, file=paste0('data/rData/', gse, '_GPL5188_Sample_Information.RDS'))


###
gse <- 'GSE54460'
seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE
phenoData <- pData(seriesMatrix[[1]])
View(phenoData)
colnames(phenoData)

phenoData <- phenoData[,c(1,2,8,42,56:70)]

colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
colnames(phenoData) <- gsub(' ', '_', colnames(phenoData))

colnames(phenoData)[3] <- 'tissue'

colnames(phenoData)[5:6] <- c('mapped_percent','bcr_status')
colnames(phenoData)[10] <- c('maybe_gleason_grade')

colnames(phenoData)[12] <- c('months_to_last_follow_up')
colnames(phenoData)[14] <- c('no_of_total_reads')
colnames(phenoData)[15] <- c('preop_psa')
colnames(phenoData)[16] <- c('pathological_t_stage')
colnames(phenoData)[19] <- c('gleason_score')



###

#GSE70770: Metadata also available in data/rData/Suppl.Table2

gse <- 'GSE70768'
seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE
phenoData <- pData(seriesMatrix[[1]])
View(phenoData)

phenoData <- phenoData[,c(1,2,8,33,38,46:59)]
colnames(phenoData)

colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
colnames(phenoData) <- gsub(' ', '_', colnames(phenoData))

colnames(phenoData)[3] <- 'tissue'
colnames(phenoData)[6:19] <- c('iclusterplus_group','age','bcr_status',
                               'clinical_stage', 'ece', 'pathological_stage','psm','preop_psa',
                               'sample_type','months_to_bcr', 'tmprss2', 'months_to_last_follow_up',
                               'tumor_percent','gleason_score')

table(phenoData$sample_type)



###
gse <- 'GSE70769'
seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE
phenoData <- pData(seriesMatrix[[1]])
View(phenoData)

phenoData <- phenoData[,c(1,2,8,34,41:51)]
colnames(phenoData)

colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
colnames(phenoData) <- gsub(' ', '_', colnames(phenoData))

colnames(phenoData)[3] <- 'tissue'
colnames(phenoData)[5:15] <- c('bcr_status','clinical_stage', 'iclusterplus_group',
                               'ece', 'pathological_stage','psm','preop_psa',
                               'months_to_bcr', 'months_to_last_follow_up',
                               'tumor_percent','gleason_score')



###
gse <- 'GSE94767'
seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE
phenoData <- pData(seriesMatrix[[1]])
View(phenoData)
colnames(phenoData)

phenoData <- phenoData[,c(1,2,8,9,32,39:45)]

colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
colnames(phenoData) <- gsub(' ', '_', colnames(phenoData))

colnames(phenoData)[3] <- 'tissue'

colnames(phenoData)[4] <- c('organism')
colnames(phenoData)[6] <- c('center')

table(phenoData$material_type)

saveRDS(phenoData, file=paste0('data/rData/', gse, '_Sample_Information.RDS'))



###
gse <- 'GSE116918'
seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE
phenoData <- pData(seriesMatrix[[1]])
View(phenoData)

phenoData <- phenoData[,c(1,2,8,29,36:43)]
colnames(phenoData)

colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
colnames(phenoData) <- gsub(' ', '_', colnames(phenoData))

colnames(phenoData)[3] <- 'tissue'
colnames(phenoData)[5:12] <- c('bcr_status','months_to_bcr',
                               'months_to_met','gleason_score', 
                               'met_status', 'age', 'preop_psa','t_stage')


###
gse <- 'GSE25136'
seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE
phenoData <- pData(seriesMatrix[[1]])
View(phenoData)
colnames(phenoData)

phenoData <- phenoData[,c(1,2,8,19,25,32)]

colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
colnames(phenoData) <- gsub(' ', '_', colnames(phenoData))

colnames(phenoData)[3] <- 'tissue'



saveRDS(phenoData, file=paste0('data/rData/', gse, '_Sample_Information.RDS'))



###
gse <- 'GSE44353'
seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE
phenoData <- pData(seriesMatrix[[1]])
View(phenoData)
colnames(phenoData)

phenoData <- phenoData[,c(1,2,8,21,27,36:38)]

colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
colnames(phenoData) <- gsub(' ', '_', colnames(phenoData))

colnames(phenoData)[3] <- 'tissue'
colnames(phenoData)

colnames(phenoData)[4] <- 'description'
colnames(phenoData)[6] <- 'bcr_status'
colnames(phenoData)[8] <- 'pathological_t_stage'

saveRDS(phenoData, file=paste0('data/rData/', gse, '_Sample_Information.RDS'))




###
gse <- 'GSE41408'
seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE
phenoData <- pData(seriesMatrix[[1]])
View(phenoData)
colnames(phenoData)

phenoData <- phenoData[,c(1,2,8,33,40:44)]

colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
colnames(phenoData) <- gsub(' ', '_', colnames(phenoData))

colnames(phenoData)[3] <- 'tissue'
colnames(phenoData)[5:9] <- c('met_status','gleason_score','preop_psa','bcr_status','pathological_t_stage')


###
gse <- 'GSE46691'

seriesMatrix <- getGEO(filename = paste0('data/fromGEO/', gse, '_series_matrix.txt.gz'))
phenoData <- pData(seriesMatrix)
View(phenoData)

phenoData <- phenoData[,c('title','geo_accession','source_name_ch1','description',
                          'contact_institute','gleason score:ch1','metastatic event:ch1')]

colnames(phenoData)
phenoData

colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
colnames(phenoData) <- gsub(' ', '_', colnames(phenoData))

colnames(phenoData)[3] <- 'tissue'



###
gse <- 'GSE51066'
seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE
phenoData <- pData(seriesMatrix[[1]])
View(phenoData)
colnames(phenoData)

phenoData <- phenoData[,c(1,2,8,21,30,38)]

colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
colnames(phenoData) <- gsub(' ', '_', colnames(phenoData))

colnames(phenoData)[3] <- 'tissue'
colnames(phenoData)

saveRDS(phenoData, file=paste0('data/rData/', gse, '_Sample_Information.RDS'))



###
gse <- 'GSE37199'
seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE
phenoData <- pData(seriesMatrix[[1]])
View(phenoData)
colnames(phenoData)

phenoData <- phenoData[,c(1,2,8,32,41:45)]

colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
colnames(phenoData) <- gsub(' ', '_', colnames(phenoData))

colnames(phenoData)[3] <- 'tissue'
colnames(phenoData)

colnames(phenoData)[5] <- 'center'
colnames(phenoData)[9] <- 'prognosis_status'

saveRDS(phenoData, file=paste0('data/rData/', gse, '_Sample_Information.RDS'))


##########################################################################################
# ===================================== Expression ===================================== #
##########################################################################################

# http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp
# R Source Package: O

# Affymetrix Human Exon 1.0 ST Array: pd.huex10st.hs.gencodeg
# Affymetrix Human Gene 2.0 ST Array: pd.huex20st.hs.gencodeg
# Affymetrix Human Transcriptome Array 2.0: pd.hta20.hs.gencodeg
# Affymetrix Human Genome U133A Array: pd.hgu133a.hs.gencodeg
# Affymetrix Human Genome U133 Plus 2.0 Array: pd.hgu133plus2.hs.gencodeg

# library(pd.hg.u133.plus.2) # from Bioconductor (NOT IN USE)

#install.packages("http://mbni.org/customcdf/24.0.0/gencodeg.download/pd.huex10st.hs.gencodeg_24.0.0.tar.gz",
#                 repos = NULL, type = "source")

library(pd.huex10st.hs.gencodeg)

gse <- 'GSE46691'

filePaths = getGEOSuppFiles(gse, baseDir = 'data/fromGEO', makeDirectory = FALSE,filter_regex = 'RAW')
untar(paste0('data/fromGEO/', gse, '_RAW.tar'), exdir = paste0('data/fromGEO/', gse, '_RAW'))

celFiles = list.celfiles(paste0('data/fromGEO/', gse, '_RAW'), full.names=T, listGzipped=T)
celFiles

rawData = read.celfiles(celFiles, pkgname = 'pd.huex10st.hs.gencodeg')

probesetData = oligo::rma(rawData)

exprData = exprs(probesetData)
colnames(exprData)
rownames(exprData)
colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
rownames(exprData) <- unlist(lapply(rownames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
colnames(exprData) <- gsub('_GBX', '', colnames(exprData))

saveRDS(exprData, file=paste0('data/rData/', gse, '_Expression_CDF24_GENCODE_RMA.RDS'))




######
gse <- 'GSE51066'

filePaths = getGEOSuppFiles(gse, baseDir = 'data/fromGEO', makeDirectory = FALSE,filter_regex = 'RAW')
untar(paste0('data/fromGEO/', gse, '_RAW.tar'), exdir = paste0('data/fromGEO/', gse, '_RAW'))

celFiles = list.celfiles(paste0('data/fromGEO/', gse, '_RAW'), full.names=T, listGzipped=T)
celFiles

rawData = read.celfiles(celFiles, pkgname = 'pd.huex10st.hs.gencodeg') #pd.huex.1.0.st.v2

probesetData = oligo::rma(rawData)

exprData = exprs(probesetData)
exprData[1:5,1:5]
colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
rownames(exprData) <- unlist(lapply(rownames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
colnames(exprData) <- gsub('_GBX', '', colnames(exprData))

saveRDS(exprData, file=paste0('data/rData/', gse, '_Expression_CDF24_GENCODE_RMA.RDS'))



########
gse <- 'GSE41408'

filePaths = getGEOSuppFiles(gse, baseDir = 'data/fromGEO', makeDirectory = FALSE,filter_regex = 'RAW')
untar(paste0('data/fromGEO/', gse, '_RAW.tar'), exdir = paste0('data/fromGEO/', gse, '_RAW'))

celFiles = list.celfiles(paste0('data/fromGEO/', gse, '_RAW'), full.names=T, listGzipped=T)
celFiles

rawData = read.celfiles(celFiles, pkgname = 'pd.huex10st.hs.gencodeg') #pd.huex.1.0.st.v2

probesetData = oligo::rma(rawData)

exprData = exprs(probesetData)
exprData[1:5,1:5]
colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
rownames(exprData) <- unlist(lapply(rownames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '_', fixed=T)[[1]][1]))

saveRDS(exprData, file=paste0('data/rData/', gse, '_Expression_CDF24_GENCODE_RMA.RDS'))




########
gse <- 'GSE116918'

seriesMatrix <- getGEO(filename = paste0('data/fromGEO/', gse, '_series_matrix.txt.gz'))
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]


annoData <- seriesMatrix@featureData@data
dim(annoData)
unique(annoData$`Ensembl Gene ID`)
View(annoData)

sum(rownames(annoData) == rownames(exprData))

filter <- which(!startsWith(annoData$`Ensembl Gene ID`, 'ENSG'))
annoData <- annoData[-filter,]
exprData <- exprData[-filter,]
dim(annoData)

filter <- grep('/', annoData$`Ensembl Gene ID`)
annoData <- annoData[-filter,]
exprData <- exprData[-filter,]
dim(annoData)

exprData <- data.frame(exprData, stringsAsFactors = F)
exprData[1:5,1:5]

exprData$ID <- annoData$`Ensembl Gene ID`
exprData$ID

exprData <- selectProbeFun(exprData)
dim(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_Expression_CDF24_GENCODE_RMA.RDS'))

#====================================================================================================================

anno <- read.table('data/fromGEO/GPL25318_ADXPCv1a520642.annot.txt', stringsAsFactors = F, header = T, sep = '\t')
anno <- read.table('data/fromGEO/GPL25318-6052.txt', stringsAsFactors = F, header = T, sep = '\t', skip = 1:100)

#BiocManager::install('pd.adxpcv1a520642') # NOT AVAILABLE

BiocManager::install('pdInfoBuilder')
library("pdInfoBuilder")
cdfFile <- "data/fromGEO/GPL25318_ADXPCv1a520642.CDF"
csvAnno <- "data/fromGEO/GPL25318_ADXPCv1a520642.annot.txt"
## csvSeq <- "Mapping250K_Nsp_probe_tab"
## 
pkg <- new("adxpcv1a520642", 
           version="0.0.1",
           author="Ruidong Li", email="bioinfo.dong@gmail.com",
           biocViews="AnnotationData",
           genomebuild="NCBI Build 37",
           cdfFile=cdfFile, 
           csvAnnoFile=csvAnno)#, 
           #csvSeqFile=csvSeq)
## 
## makePdInfoPackage(pkg, destDir=".")


filePaths = getGEOSuppFiles(gse, baseDir = 'data/fromGEO', makeDirectory = FALSE,filter_regex = 'RAW')
untar(paste0('data/fromGEO/', gse, '_RAW.tar'), exdir = paste0('data/fromGEO/', gse, '_RAW'))

celFiles = list.celfiles(paste0('data/fromGEO/', gse, '_RAW'), full.names=T, listGzipped=T)
celFiles

rawData = read.celfiles(celFiles) #pd.huex.1.0.st.v2

probesetData = oligo::rma(rawData)

exprData = exprs(probesetData)
exprData[1:5,1:5]
colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
rownames(exprData) <- unlist(lapply(rownames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '_', fixed=T)[[1]][1]))

saveRDS(exprData, file=paste0('data/rData/', gse, '_Expression_CDF24_GENCODE_RMA.RDS'))
#===============================================================================================================



########
gse <- 'GSE21034'

filePaths = getGEOSuppFiles(gse, baseDir = 'data/fromGEO', makeDirectory = FALSE,filter_regex = 'RAW')
untar(paste0('data/fromGEO/', gse, '_RAW.tar'), exdir = paste0('data/fromGEO/', gse, '_RAW'))

celFiles = list.celfiles(paste0('data/fromGEO/', gse, '_RAW'), full.names=T, listGzipped=T)
celFiles

celFiles <- celFiles[grep('GSM526', celFiles)]

rawData <- read.celfiles(celFiles, pkgname = 'pd.huex10st.hs.gencodeg') #pd.huex.1.0.st.v2

probesetData <- oligo::rma(rawData)

exprData <- exprs(probesetData)
exprData[1:5,1:5]
#colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
rownames(exprData) <- unlist(lapply(rownames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '_', fixed=T)[[1]][1]))

saveRDS(exprData, file=paste0('data/rData/', gse, '_Expression_CDF24_GENCODE_RMA.RDS'))



###################

## GSE54460

# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz

# STAR, featureConts

countMatrix <- read.table('data/fromSRA/GSE54460/GSE54460_featureCounts_Gencode.v32.txt',
                       header = T, stringsAsFactors = F)
countMatrix

countMatrix$Geneid <- unlist(lapply(countMatrix$Geneid, function(x) strsplit(x, '.', fixed=T)[[1]][1]))
filter <- which(duplicated(countMatrix$Geneid))

countMatrix <- countMatrix[-filter,]
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[,-c(1:2)]
colnames(countMatrix) <- gsub('BAM.|.bam', '', colnames(countMatrix))
dim(countMatrix)

dge <-  DGEList(counts = countMatrix)

### TMM normalization
dge = calcNormFactors(dge, method = 'TMM')

exprLogCPM <- edgeR::cpm(dge,log = TRUE) ### for visualization
dim(exprLogCPM)

saveRDS(exprLogCPM, file='data/rData/GSE54460_Expression_LogCPM_All_Genes.RDS')






###################

### GPC-GENE, 2017; GSE107299

gse <- 'GSE107299'

library(pd.hugene20st.hs.gencodeg)
library(pd.hta20.hs.gencodeg)

filePaths = getGEOSuppFiles(gse, baseDir = 'data/fromGEO', makeDirectory = FALSE,filter_regex = 'RAW')
untar(paste0('data/fromGEO/', gse, '_RAW.tar'), exdir = paste0('data/fromGEO/', gse, '_RAW'))

celFiles = list.celfiles(paste0('data/fromGEO/', gse, '_RAW'), full.names=T, listGzipped=T)
celFiles

batch2 <- readRDS(file=paste0('data/rData/', gse, '_GPL17586_Sample_Information.RDS'))
batch2
dim(batch2)

# BATCH 2, HTA2.0, GPL17586
idx <- which(unlist(lapply(strsplit(celFiles, '_|/', fixed=F), function(x) x[5])) %in% 
               paste0('GSM', 2981313:2981378))
idx
celFiles <- celFiles[idx]
rawData <- read.celfiles(celFiles, pkgname = 'pd.hta20.hs.gencodeg') #pd.huex.1.0.st.v2

# OTHER BATCHES, HuG2.0
idx <- which(!unlist(lapply(strsplit(celFiles, '_|/', fixed=F), function(x) x[5])) %in% 
               paste0('GSM', 2981313:2981378))
celFiles <- celFiles[idx]
celFiles
rawData <- read.celfiles(celFiles, pkgname = 'pd.hugene20st.hs.gencodeg') #pd.huex.1.0.st.v2


probesetData <- oligo::rma(rawData)

exprData <- exprs(probesetData)
exprData[1:5,1:5]
dim(exprData)
exprData <- exprData[grep('ENSG', rownames(exprData)),]

View(exprData)
#colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
rownames(exprData) <- unlist(lapply(rownames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '_', fixed=T)[[1]][1]))

saveRDS(exprData, file=paste0('data/rData/', gse, '_GPL17586_Expression_CDF24_GENCODE_RMA.RDS'))
saveRDS(exprData, file=paste0('data/rData/', gse, '_GPL16686_Expression_CDF24_GENCODE_RMA.RDS'))



##################
### GPC-GENE, 2017; GSE107299
### From GEO
gse <- 'GSE107299'

# ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
# ftp://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz

ensembl <- readRDS(file='data/ENSEMBL_Annotation_Human_V98.RDS')
gencode <- readRDS(file='data/GENCODE_Annotation_Human_V32.RDS')
hgnc <- read.table(file='data/hgnc_complete_set.txt', sep='\t', quote='', comment.char='', 
                   header = T, stringsAsFactors = F)

ncbi.gene.info <- read.table('data/Annotation/NCBI_Homo_sapiens.gene_info.20191230.gz', header = T, stringsAsFactors = F, 
                             sep='\t', quote='', comment.char = '')

ncbi.gene.info$Ensembl <- str_extract(ncbi.gene.info$dbXrefs, 'ENSG\\d+')

ncbi.gene.history <- read.table(file='data/Annotation/NCBI_Homo_sapiens.gene_history.20191230.txt', sep='\t', header=T,
                                stringsAsFactors = F, quote = '', comment.char = '')

idx<- which(is.na(ncbi.gene.info$Ensembl))
ncbi.gene.info$Ensembl[idx] <- hgnc$ensembl_gene_id[match(ncbi.gene.info$GeneID[idx], hgnc$entrez_id)]

idx<- which(is.na(ncbi.gene.info$Ensembl))
ncbi.gene.info$Ensembl[idx] <- ensembl$ensembl_gene_id[match(ncbi.gene.info$GeneID[idx], ensembl$entrezgene_id)]

idx<- which(is.na(ncbi.gene.info$Ensembl))
idx

ncbi.gene.info[idx,]

gtf <- getGENCODEAnnotation(gtf.file = 'data/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz')
gtf[1:5,]

gtf <- readGFF('data/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz')
gtf[1:5,]


###
exprData <- read.table(file='data/fromGEO/GSE107299_Matrix_processed_data.tsv',stringsAsFactors = F,
                       header = T, sep = '\t', quote='')
exprData

idx <- which(exprData$GeneID %in% ncbi.gene.history$Discontinued_GeneID)
idx
exprData$GeneID[idx] <- ncbi.gene.history$GeneID[match(exprData$GeneID[idx], ncbi.gene.history$Discontinued_GeneID)]

exprData$GeneID[!exprData$GeneID %in% ncbi.gene.info$GeneID]

exprData <- exprData[-which(exprData$GeneID=='-'),]
rownames(exprData) <- 1:nrow(exprData)
View(exprData)

dupids <- exprData$GeneID[duplicated(exprData$GeneID)]
dupids

View(exprData[which(exprData$GeneID %in% dupids),])

filter <- c(13332,1264,3979,344,18880,291,22088,18805,13461,19254,
            21106,1299,3959,3977,2184,6779,10190,23402,3882,9693,
            8514,546,10506,3193,1319,3908,3978,3975,153,3941,3447,
            13896,3214,2754,13305,12776,3730,23783,10363,2450,3906,
            6699,2423,18719,1481,7465,3916,1308,12163,3071,1223,13423,
            18648,3396,3971,1518,9499,18227,3961)

exprData <- exprData[-filter,]

exprData <- add_column(.data = exprData, .after = 1, Ensembl=NA)
exprData[1:5,1:5]

dim(exprData)
sum(exprData$GeneID %in% ensembl$entrezgene_id)
sum(exprData$GeneID %in% hgnc$entrez_id)
sum(exprData$GeneID %in% ncbi.gene.info$GeneID)

idx <- match(exprData$GeneID, ncbi.gene.info$GeneID)
idx

exprData$Ensembl <- ncbi.gene.info$Ensembl[idx]



final.anno <- readRDS('data/Annotation/Homo_Sapiens_Gene_Annotation_ENSEMBL_HGNC_ENTREZ.RDS')
idx <- match(exprData$GeneID, final.anno$entrez_id)
idx


exprData <- add_column(.data = exprData, .after = 1, Name=NA)
exprData[1:5,1:5]

exprData$Name <- ifelse(is.na(ncbi.gene.info[idx,]$Ensembl) | ncbi.gene.info[idx,]$Ensembl=='',
                        ncbi.gene.info[idx,]$Symbol, ncbi.gene.info[idx,]$Ensembl)

which(is.na(final.anno[idx,]$ensembl_id))


View(exprData)

rownames(exprData) <- exprData$Name
exprData <- exprData[,-c(1:9)]
colnames(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_Merge_Processed_Expression_From_GEO.RDS'))

###############################################################

### GSE70700: GSE70768, GSE70769

# GSE70768

#=========================================================================================================
# from ENSEMBL
attributes <- c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name', 'illumina_humanht_12_v4')
ILLUMINA.HumanHT.12.V4 <- getENSEMBLAnnotation(attributes = attributes)

View(ILLUMINA.HumanHT.12.V4)
dim(ILLUMINA.HumanHT.12.V4)

sum(ILLUMINA.HumanHT.12.V4$entrezgene_id %in% ncbi.gene.history$Discontinued_GeneID) #51

View(ILLUMINA.HumanHT.12.V4[ILLUMINA.HumanHT.12.V4$entrezgene_id %in% ncbi.gene.history$Discontinued_GeneID,])

sum(startsWith(ILLUMINA.HumanHT.12.V4$illumina_humanht_12_v4,'ILMN'))
idx <- which(startsWith(ILLUMINA.HumanHT.12.V4$illumina_humanht_12_v4,'ILMN'))
unique(ILLUMINA.HumanHT.12.V4$ensembl_gene_id[idx])
length(unique(ILLUMINA.HumanHT.12.V4$illumina_humanht_12_v4))
# 38943 probes were annotated to 30254 ensembl genes
# one probe may be annotated to multiple genes, e.g., ILMN_3311070

idx <- which(is.na(ILLUMINA.HumanHT.12.V4$entrezgene_id))
idx

ids <- ILLUMINA.HumanHT.12.V4$ensembl_gene_id[idx]
ids

ILLUMINA.HumanHT.12.V4$entrezgene_id[idx] <- ncbi.gene.info$GeneID[match(ids, ncbi.gene.info$Ensembl)]


idx <- which(ILLUMINA.HumanHT.12.V4$entrezgene_id %in% ncbi.gene.history$Discontinued_GeneID)
idx

ids <- ILLUMINA.HumanHT.12.V4$entrezgene_id[idx]
ids
ids <- ncbi.gene.history$GeneID[match(ids, ncbi.gene.history$Discontinued_GeneID)]

ILLUMINA.HumanHT.12.V4$entrezgene_id[idx] <- ifelse(ids=='-', NA, ids)

View(ILLUMINA.HumanHT.12.V4)



t <- ILLUMINA.HumanHT.12.V4 %>% group_by(illumina_humanht_12_v4) %>%
  summarise(ensembl=paste(ensembl_gene_id, collapse=';'))
View(t)
colnames(ILLUMINA.HumanHT.12.V4)

length(unique(t$ensembl))

grep(';',t$ensembl)


idx <- which(ILLUMINA.HumanHT.12.V4$ensembl_gene_id %in% gencode$ensembl)
ILLUMINA.HumanHT.12.V4 <- ILLUMINA.HumanHT.12.V4[idx,]



# from Bioconductor
BiocManager::install("illuminaHumanv4.db")
library(illuminaHumanv4.db)

help('select')
columns(illuminaHumanv4.db)

x <- select(illuminaHumanv4.db, 
            keys = rownames(exprData), 
            columns=c("SYMBOL","ENTREZID",'ENSEMBL'), 
            keytype="PROBEID")
View(x)

sum(!is.na((x$ENSEMBL)))
sum(!is.na((x$ENTREZID)))
sum(!is.na((x$SYMBOL)))
unique(x$ENTREZID)
unique(x$ENSEMBL)
sum(x$ENTREZID %in% ncbi.gene.history$Discontinued_GeneID) # 25
dim(x)
sum(x$ENTREZID %in% ncbi.gene.info$GeneID)

dim(exprData)
sum(rownames(exprData) %in% x$PROBEID)


length(unique(x$PROBEID[which(!is.na((x$ENTREZID)))])) 
# 35722 probes were annotated to 21631 entrez genes, and 24038 ensembl genes



# from GEO (IN USE)
annoData <- read.table('data/fromGEO/GPL10558-50081.txt', header = T, sep = '\t', 
                       comment.char = '', quote='', skip = 30, stringsAsFactors = F)

which(is.na(annoData$Entrez_Gene_ID))

idx <- which(annoData$Entrez_Gene_ID %in% ncbi.gene.history$Discontinued_GeneID)
idx

ids <- annoData$Entrez_Gene_ID[idx]
ids

annoData$Entrez_Gene_ID[idx] <- ncbi.gene.history$GeneID[match(ids, ncbi.gene.history$Discontinued_GeneID)]
ncbi.gene.history$GeneID[match(ids, ncbi.gene.history$Discontinued_GeneID)]

filter <- which(annoData$Entrez_Gene_ID=='-' | 
                  is.na(annoData$Entrez_Gene_ID) | 
                  annoData$Entrez_Gene_ID %in% ncbi.gene.history$Discontinued_GeneID)
filter
annoData <- annoData[-filter, ]

dim(annoData)
View(annoData)

annoData <- add_column(.data = annoData, .after = 10, Ensembl=NA)

idx <- match(annoData$Entrez_Gene_ID, ncbi.gene.info$GeneID)
idx

annoData$Ensembl <- ncbi.gene.info$Ensembl[idx]

unique(annoData$Ensembl)
unique(annoData$Entrez_Gene_ID)

View(annoData[which(!annoData$Ensembl %in% t$ensembl),])
View(t[which(!t$ensembl %in% annoData$Ensembl),])

# 37498 probes were annotated to entrez genes


keep <- which(startsWith(annoData$Ensembl, 'ENSG'))
keep

annoData <- annoData[keep,]
rownames(annoData) <- annoData$ID

#=========================================================================================================

gse <- 'GSE70768'
gse <- 'GSE70769'
seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE

exprData <- exprs(seriesMatrix[[1]])
exprData <- data.frame(exprData[annoData$ID,])
exprData$ID <- annoData$Ensembl
View(exprData)

idx <- which(!is.na(rowSums(exprData[,-ncol(exprData)])))
idx
exprData[idx,]

exprData <- selectProbeFun(exprData[idx,])

saveRDS(exprData, file=paste0('data/rData/', gse, '_Expression_from_GEO.RDS'))


#######################

gse <- 'GSE25136'

library(pd.hgu133a.hs.gencodeg)

filePaths = getGEOSuppFiles(gse, baseDir = 'data/fromGEO', makeDirectory = FALSE,filter_regex = 'RAW')
untar(paste0('data/fromGEO/', gse, '_RAW.tar'), exdir = paste0('data/fromGEO/', gse, '_RAW'))

celFiles = list.celfiles(paste0('data/fromGEO/', gse, '_RAW'), full.names=T, listGzipped=T)
celFiles

rawData = read.celfiles(celFiles, pkgname = 'pd.hgu133a.hs.gencodeg')

probesetData = oligo::rma(rawData)

exprData = exprs(probesetData)
exprData[1:5,1:5]
colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
rownames(exprData) <- unlist(lapply(rownames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
#colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '_', fixed=T)[[1]][1]))

exprData <- exprData[which(startsWith(rownames(exprData), 'ENSG')),]
dim(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_Expression_CDF24_GENCODE_RMA.RDS'))


#######################

gse <- 'GSE94767'

library(pd.huex10st.hs.ensg)

filePaths = getGEOSuppFiles(gse, baseDir = 'data/fromGEO', makeDirectory = FALSE, filter_regex = 'RAW')
untar(paste0('data/fromGEO/', gse, '_RAW.tar'), exdir = paste0('data/fromGEO/', gse, '_RAW'))

celFiles = list.celfiles(paste0('data/fromGEO/', gse, '_RAW'), full.names=T, listGzipped=T)
celFiles <- celFiles[1:236]

rawData = read.celfiles(celFiles, pkgname = 'pd.huex10st.hs.ensg')

probesetData = oligo::rma(rawData)

exprData = exprs(probesetData)
exprData[1:5,1:5]
#colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
rownames(exprData) <- unlist(lapply(rownames(exprData), function(x) strsplit(x, '_', fixed=T)[[1]][1]))
colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '_', fixed=T)[[1]][1]))

which(!startsWith(rownames(exprData), 'ENSG'))
exprData <- exprData[which(startsWith(rownames(exprData), 'ENSG')),]
dim(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_Expression_CDF24_GENCODE_RMA.RDS'))

##############

# E-MTAB-6128

gse <- 'E_MTAB_6128'

for (i in 1:8) {
  system(paste0('wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6128/E-MTAB-6128.raw.', i, 
                '.zip -P data/fromArrayExpress/E-MTAB-6128/'))
  unzip(paste0('data/fromArrayExpress/E-MTAB-6128/E-MTAB-6128.raw.', i, '.zip'), exdir = 'data/fromArrayExpress/E-MTAB-6128/')
}


library(pd.hugene20st.hs.gencodeg)

celFiles = list.celfiles(paste0('data/fromArrayExpress/E-MTAB-6128/'), full.names=T, listGzipped=T)
celFiles

rawData = read.celfiles(celFiles, pkgname = 'pd.hugene20st.hs.gencodeg')

probesetData = oligo::rma(rawData)

exprData = exprs(probesetData)
exprData[1:5,1:5]
colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
rownames(exprData) <- unlist(lapply(rownames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
#colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '_', fixed=T)[[1]][1]))

which(!startsWith(rownames(exprData), 'ENSG'))
exprData <- exprData[which(startsWith(rownames(exprData), 'ENSG')),]
dim(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_Expression_CDF24_GENCODE_RMA.RDS'))



##################

###
gse <- 'GSE44353'
seriesMatrix <- getGEO(gse, AnnotGPL = TRUE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE

exprData <- exprs(seriesMatrix[[1]])
View(exprData)
annoData <- seriesMatrix[[1]]@featureData@data
View(annoData)

idx <- which(annoData$Entrez_Gene_ID %in% ncbi.gene.history$Discontinued_GeneID)
ids <- annoData$Entrez_Gene_ID[idx]

ncbi.gene.history$GeneID[match(ids, ncbi.gene.history$Discontinued_GeneID)]

annoData <- add_column(.data = annoData, .after = 5, Ensembl=NA)

idx <- match(annoData$Entrez_Gene_ID, ncbi.gene.info$GeneID)
annoData$Ensembl <- ncbi.gene.info$Ensembl[idx]

View(annoData)

annoData <- annoData[-which(is.na(annoData$Ensembl)),]

exprData <- exprData[rownames(annoData),]
dim(exprData)

annoData[which(duplicated(annoData$Ensembl)),]

exprData <- data.frame(exprData)
exprData$ID <- annoData$Ensembl

exprData <- selectProbeFun(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_Expression_from_GEO.RDS'))


#######################

## GSE37199
gse <- 'GSE37199'

library(pd.hgu133plus2.hs.gencodeg)

filePaths = getGEOSuppFiles(gse, baseDir = 'data/fromGEO', makeDirectory = FALSE, filter_regex = 'RAW')
untar(paste0('data/fromGEO/', gse, '_RAW.tar'), exdir = paste0('data/fromGEO/', gse, '_RAW'))

celFiles = list.celfiles(paste0('data/fromGEO/', gse, '_RAW'), full.names=T, listGzipped=T)
celFiles

rawData = read.celfiles(celFiles, pkgname = 'pd.hgu133plus2.hs.gencodeg')

probesetData = oligo::rma(rawData)

exprData = exprs(probesetData)
exprData[1:5,1:5]
View(exprData)
#colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
rownames(exprData) <- unlist(lapply(rownames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '_', fixed=T)[[1]][1]))

which(!startsWith(rownames(exprData), 'ENSG'))
exprData <- exprData[which(startsWith(rownames(exprData), 'ENSG')),]
dim(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_Expression_CDF24_GENCODE_RMA.RDS'))

##############################

### Metastatic Prostate Adenocarcinoma (SU2C/PCF Dream Team, PNAS 2019)

exprDataPolyA <- read.table('data/cBioPortal/prad_su2c_2019_data_mRNA_seq_fpkm_polya.txt', header=T,
                            stringsAsFactors = F)
dim(exprDataPolyA)
exprDataPolyA$Hugo_Symbol[which(duplicated(exprDataPolyA$Hugo_Symbol))]
View(exprDataPolyA)


exprDataCapture <- read.table('data/cBioPortal/prad_su2c_2019_data_mRNA_seq_fpkm_capture.txt', header=T,
                            stringsAsFactors = F, row.names = 1)
dim(exprDataCapture)

exprDataCapture$Hugo_Symbol[which(duplicated(exprDataCapture$Hugo_Symbol))]
exprDataCapture[1:5,1:5]

intersect(colnames(exprDataPolyA), colnames(exprDataCapture))
View(exprDataCapture)

exprDataPolyA$Hugo_Symbol==exprDataCapture$Hugo_Symbol


###########################################################################################################






########################################################################################
# ==================================== ExpressionSet ==================================#
########################################################################################

gse <- 'GSE46691'
gse <- 'GSE51066'
gse <- 'GSE41408'
gse <- 'GSE59745'
gse <- 'GSE116918'
gse <- 'GSE25136'
gse <- 'GSE37199'


exprData <- readRDS(file=paste0('data/rData/', gse, '_Expression_CDF24_GENCODE_RMA.RDS'))
phenoData <- readRDS(file=paste0('data/rData/', gse, '_Sample_Information.RDS'))
View(phenoData)

gse <- 'GSE70768'
gse <- 'GSE70769'

exprData <- readRDS(file=paste0('data/rData/', gse, '_Expression_from_GEO.RDS'))
phenoData <- readRDS(file=paste0('data/rData/', gse, '_Sample_Information.RDS'))
View(phenoData)

rownames(phenoData) == colnames(exprData)

#phenoData <- read_xlsx('data/rData/Ross_Adams_2015_Suppl_Table2_Associated clinical metadata for Cambridge and Stockholm discovery and validation cohorts.xlsx')
#phenoData

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/rData/', gse, '_CDF24_GENCODE_RMA_eSet.RDS'))

saveRDS(eSet, file=paste0('data/rData/', gse, '_Expression_from_GEO_eSet.RDS'))

View(phenoData)



#############
gse <- 'GSE21034'

exprData <- readRDS(file=paste0('data/rData/', gse, '_Expression_CDF24_GENCODE_RMA.RDS'))

phenoDataGEO <- readRDS(file=paste0('data/rData/', gse, '_GPL5188_Sample_Information.RDS'))
dim(phenoDataGEO)

rownames(phenoDataGEO) == colnames(exprData)

View(phenoDataGEO)

phenoData <- read.table('data/cBioPortal/prad_mskcc_clinical_data.tsv', stringsAsFactors = F,
                        header=T, sep='\t')

phenoData <- phenoData[match(phenoDataGEO$sample_id, phenoData$Sample.ID),]

phenoData$Sample.ID==phenoDataGEO$sample_id

rownames(phenoData) <- phenoDataGEO$geo_accession

rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/rData/', gse, '_CDF24_GENCODE_RMA_eSet.RDS'))



##########################

dataset <- 'DKFZ2018'

exprData <- read.table('data/cBioPortal/prostate_dkfz_2018_RNA_Seq_expression_median.txt', sep='\t', 
                       header = T, stringsAsFactors = F)

phenoData <- read.table('data/cBioPortal/prostate_dkfz_2018_clinical_data.tsv', stringsAsFactors = F,
                        header=T, sep='\t')

rownames(phenoData) <- phenoData$Sample.ID

exprData[1:5,1:5]
dim(exprData)

## ENSEMBL 62; ftp://ftp.ensembl.org/pub/release-62/gtf/homo_sapiens/
gtf <- readGFF('data/Annotation/Homo_sapiens.GRCh37.62.gtf.gz', version=2L)
filter <- which(duplicated(gtf$gene_id))
gtf <- gtf[-filter,]
gtf

exprData <- add_column(.data = exprData, .before = 3, Ensembl=NA)
exprData$Ensembl <- gtf$gene_id[match(exprData$Hugo_Symbol, gtf$gene_name)]
exprData[1:5,1:5]

filter <- which(duplicated(exprData$Ensembl))
filter

exprData <- exprData[-filter,]
rownames(exprData) <-exprData$Ensembl

exprData <- exprData[,-c(1:3)]

ovlp <- intersect(rownames(phenoData), colnames(exprData))
ovlp

exprData <- exprData[,ovlp]
phenoData <- phenoData[ovlp,]

rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/rData/', dataset, '_cBioPortal_RPKM_eSet.RDS'))

###################

### GPC-GENE, 2017; GSE107299

gse <- 'GSE107299'

#================================================================================
phenoDataGEO1 <- readRDS(file=paste0('data/rData/', gse, '_GPL17586_Sample_Information.RDS'))
phenoDataGEO2 <- readRDS(file=paste0('data/rData/', gse, '_GPL16686_Sample_Information.RDS'))

phenoDataGEO <- rbind(phenoDataGEO1, phenoDataGEO2)
dim(phenoDataGEO)

rownames(phenoDataGEO) <- gsub('prostate tumour ', '', phenoDataGEO$title)

rownames(phenoDataGEO) == colnames(exprData)
phenoDataGEO <- phenoDataGEO[colnames(exprData),]


View(phenoDataGEO)

phenoDataGEO$sample_id <- gsub('prostate tumour ', '', phenoDataGEO$title)

phenoDataFromcBio <- read.table('data/cBioPortal/prad_cpcg_2017_clinical_data.tsv', stringsAsFactors = F,
                        header=T, sep='\t', quote='', comment.char = '')
phenoDataFromcBio
table(phenoDataFromcBio$Cohort)

View(phenoDataFromcBio)

phenoDataFromcBio$Sample.ID <- gsub('-F1', '', phenoDataFromcBio$Sample.ID)

phenoDataFromcBio <- phenoDataFromcBio[,c('Patient.ID','Sample.ID','Diagnosis.Age',)]
colnames(phenoDataFromcBio)


phenoDataFromPaper <- read.table(file = 'data/fromGEO/GSE107299_Clinical_Data_From_Paper.txt', header=T, 
                        sep='\t', stringsAsFactors = F)

View(phenoDataFromPaper)


samples <- unique(c(phenoDataFromcBio$Sample.ID, phenoDataFromPaper$Patient_ID))
samples

phenoDataFromcBio <- phenoDataFromcBio[match(samples, phenoDataFromcBio$Sample.ID), ]
phenoDataFromPaper <- phenoDataFromPaper[match(samples, phenoDataFromPaper$Patient_ID), ]

phenoData <- cbind(phenoDataFromcBio, phenoDataFromPaper)
View(phenoData)

colnames(phenoData)

phenoData[which(phenoData$Patient.ID=='CPCG0102-F1'),]

colnames(phenoData)

phenoData <- phenoData[,c(1,3,4,5, 10,11,12,15,19,20,21,24,25,28,32,33,34,58,59,63,70,71,75,77,82,83,91,
                         101,102,103,107,108)]

idx <- which(is.na(phenoData$Sample.ID))
idx

phenoData$Study.ID[idx] <- 'cancer_cell_2019'
phenoData$Sample.ID[idx] <- phenoData$Patient_ID[idx]
phenoData$Diagnosis.Age[idx] <- phenoData$Age_at_Treatment[idx]

idx <- which(is.na(phenoData$Neoplasm.American.Joint.Committee.on.Cancer.Clinical.Primary.Tumor.T.Stage))
phenoData$Neoplasm.American.Joint.Committee.on.Cancer.Clinical.Primary.Tumor.T.Stage[idx] <- phenoData$Clinical_T.category[idx]


idx <- is.na(phenoData$Biochemical_Recurrence)
phenoData$Biochemical_Recurrence[idx] <- phenoData$Biochemical.Recurrence.Indicator[idx]

phenoData$Time_to_Biochemical_Recurrence_or_Max_Follow_Up <- phenoData$Time_to_Biochemical_Recurrence_or_Max_Follow_Up*365.25

idx <- is.na(phenoData$Time_to_Biochemical_Recurrence_or_Max_Follow_Up)
phenoData$Time_to_Biochemical_Recurrence_or_Max_Follow_Up[idx] <- phenoData$Days.to.biochemical.recurrence.first[idx]


idx <- which(phenoData$Biochemical_Recurrence=='NO' & phenoData$Days.to.Last.Followup>phenoData$Days.to.biochemical.recurrence.first)
phenoData$Time_to_Biochemical_Recurrence_or_Max_Follow_Up[idx] <- phenoData$Days.to.Last.Followup[idx]

idx <- which(is.na(phenoData$Time_to_Biochemical_Recurrence_or_Max_Follow_Up))
phenoData$Time_to_Biochemical_Recurrence_or_Max_Follow_Up[idx] <- phenoData$Days.to.Last.Followup[idx]

phenoData$Time_to_Biochemical_Recurrence_or_Max_Follow_Up <- round(phenoData$Time_to_Biochemical_Recurrence_or_Max_Follow_Up)
round(phenoData$Time_to_Biochemical_Recurrence_or_Max_Follow_Up)

phenoData$Biochemical_Recurrence[phenoData$Biochemical_Recurrence=='TRUE'] <- 'YES'
phenoData$Biochemical_Recurrence[phenoData$Biochemical_Recurrence=='FALSE'] <- 'NO'


View(data.frame(phenoData$Days.to.biochemical.recurrence.first, phenoData$Days.to.Last.Followup,
                phenoData$Biochemical.Recurrence.Indicator, 
                phenoData$Biochemical_Recurrence,
                phenoData$Time_to_Biochemical_Recurrence_or_Max_Follow_Up))


View(phenoData)
colnames(phenoData)
phenoData <- phenoData[,c(1,2,3,5,6,14,21,24,29,15:19,12,13,31,32,22,25)]

View(phenoData)

colnames(phenoData)[3] <- 'Diagnosis.Age'
colnames(phenoData)[4:5] <- c('Clinical.M.Stage','Clinical.T.Stage')

phenoData <- phenoData[match(phenoDataGEO$sample_id, phenoData$Sample.ID),]
phenoData

phenoData <- cbind(phenoDataGEO, phenoData)
phenoData

saveRDS(phenoData, file=paste0('data/rData/', gse, '_Organized_Meatadata.RDS'))
phenoData <- readRDS(file=paste0('data/rData/', gse, '_Organized_Meatadata.RDS'))
View(phenoData)
#================================================================================

exprData <- readRDS(file=paste0('data/rData/', gse, '_GPL17586_Expression_CDF24_GENCODE_RMA.RDS'))
exprData <- readRDS(file=paste0('data/rData/', gse, '_GPL16686_Expression_CDF24_GENCODE_RMA.RDS'))
exprData <- readRDS(file=paste0('data/rData/', gse, '_Merge_Processed_Expression_From_GEO.RDS'))
colnames(exprData)

phenoData <- readRDS(file=paste0('data/rData/', gse, '_Organized_Meatadata.RDS'))
rownames(phenoData) <- phenoData$geo_accession

rownames(phenoData) == colnames(exprData)

phenoData <- phenoData[colnames(exprData),]

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/rData/', gse, '_GPL17586_CDF24_GENCODE_RMA_eSet.RDS'))
saveRDS(eSet, file=paste0('data/rData/', gse, '_GPL16686_CDF24_GENCODE_RMA_eSet.RDS'))

saveRDS(eSet, file=paste0('data/rData/', gse, '_Merge_Prcessed_From_GEO_eSet.RDS'))

sum(!is.na(phenoData$Time_to_Biochemical_Recurrence_or_Max_Follow_Up))
sum(!is.na(phenoData$Time_to_Biochemical_Recurrence_or_Max_Follow_Up))


########

gse <- 'GSE94767'

exprData <- readRDS(file=paste0('data/rData/', gse, '_Expression_CDF24_GENCODE_RMA.RDS'))
phenoDataFromGEO <- readRDS(file=paste0('data/rData/', gse, '_Sample_Information.RDS'))
phenoDataFromGEO <- phenoDataFromGEO[-which(phenoDataFromGEO$organism=='Mus musculus'),]
View(phenoDataFromGEO)
rownames(phenoDataFromGEO) <- gsub('Prostate - ', '', phenoDataFromGEO$title)

phenoDataFromPaper <- read_xlsx('data/fromGEO/GSE94767_Clinical_Data_From_Paper_Supp.xlsx')
phenoDataFromPaper
View(phenoDataFromPaper)
phenoDataFromPaper <- data.frame(phenoDataFromPaper,stringsAsFactors = F)
rownames(phenoDataFromPaper) <- phenoDataFromPaper$Sample.ID


phenoDataFromPaper <- phenoDataFromPaper[rownames(phenoDataFromGEO),]
phenoData <- cbind(phenoDataFromGEO, phenoDataFromPaper)
rownames(phenoData) <- phenoData$geo_accession
View(phenoData)

View(exprData)
rownames(phenoData) == colnames(exprData)

#phenoData <- read_xlsx('data/rData/Ross_Adams_2015_Suppl_Table2_Associated clinical metadata for Cambridge and Stockholm discovery and validation cohorts.xlsx')
#phenoData

phenoData <- phenoData[,-11]

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/rData/', gse, '_CDF24_GENCODE_RMA_eSet.RDS'))

###############

# E-MTAB-6128
gse <- 'E_MTAB_6128'

exprData <- readRDS(file=paste0('data/rData/', gse, '_Expression_CDF24_GENCODE_RMA.RDS'))

phenoData <- read.table('data/fromArrayExpress/E-MTAB-6128/E-MTAB-6128.sdrf.txt', header = T, stringsAsFactors = F,
                        sep = '\t')

rownames(phenoData) <- phenoData$Assay.Name
phenoData <- phenoData[colnames(exprData),]

rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/rData/', gse, '_CDF24_GENCODE_RMA_eSet.RDS'))

#######
# GSE44353

gse <- 'GSE44353'

exprData <- readRDS(exprData, file=paste0('data/rData/', gse, '_Expression_from_GEO.RDS'))

phenoData <- readRDS(file=paste0('data/rData/', gse, '_Sample_Information.RDS'))
View(phenoData)

rownames(phenoData) == colnames(exprData)

#phenoData <- read_xlsx('data/rData/Ross_Adams_2015_Suppl_Table2_Associated clinical metadata for Cambridge and Stockholm discovery and validation cohorts.xlsx')
#phenoData

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/rData/', gse, '_from_GEO_eSet.RDS'))


##################
## GSE54460

gse <- 'GSE54460'

phenoDataFromPaper <- read_xlsx(path = 'data/fromGEO/GSE54460_Clinical_Data_From_Paper_Table_S2.xlsx')
#phenoDataFromPaper <- data.frame(phenoDataFromPaper, stringsAsFactors = F)
#colnames(phenoDataFromPaper) <- phenoDataFromPaper[2,]
#phenoDataFromPaper <- phenoDataFromPaper[-c(1:2),]

phenoDataFromGEO <- readRDS(file=paste0('data/rData/', gse, '_Sample_Information.RDS'))
phenoDataFromSRA <- read.table(file='data/fromSRA/GSE54460/SraRunTable.txt', sep=',', header=T,
                               stringsAsFactors = F)

idx <- match(phenoDataFromGEO$geo_accession, phenoDataFromSRA$GEO_Accession)
idx
rownames(phenoDataFromGEO) <- phenoDataFromSRA$Run[idx]
phenoData <- phenoDataFromGEO

exprData <- readRDS('data/rData/GSE54460_Expression_LogCPM_All_Genes.RDS')

rownames(phenoData) == colnames(exprData)


eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/rData/GSE54460_LogCPM_All_Genes_eSet.RDS'))





########################################################################################
# ========================== Final ExpressionSet for Database =========================#
########################################################################################

### Datasets

# TCGA-PRAD     == Cell 2015                == GDC         == RNAseq: CPM              == Time to BCR
# CPC-GENE      == Nature_2017              == GSE107299   == HuGene 2.0 ST & HTA 2.0  == Time to BCR (Download merged data from GEO)
# MSKCC         ==_Cancer_Cell_2010         == GSE21034    == HuEx 1.0 ST              == Time to BCR
# DKFZ          == Cancer Cell 2018         == cBioPortal  == RNAseq: RPKM             == Time to BCR (Download RPKM from cBioportal, one patient may have multiple samples)
# GSE54460      == Cancer Research 2014     == GSE54460    == RNAseq: CPM              == Time to BCR (Realignment)
# CamCap        == Nature Genetics 2016     == GSE70768    == Illumina HumanHT-12 V4.0 == Time to BCR
# Stockholm     == Nature Genetics 2016     == GSE70769    == Illumina HumanHT-12 V4.0 == Time to BCR
# DESNT         == Eur Urol Focus 2017      == GSE94767    == HuEx 1.0 ST              == Time to BCR
# E-MTAB-6128   == Annuals of Oncology 2018 == E-MTAB-6128 == HuGene 2.0 ST            == Time to BCR
# Belfast       == Journal of Clinical Oncology 2017  == GSE116918 == ADXPCv1a520642   == Time to BCR/Metastasis

# GSE25136      == Prostate 2019            == GSE25136    == HG-U133A                 == BCR Status
# GSE44353      == PNAS 2013 == GSE44353    == Prostate Cancer DASL Panel 1.5K         == BCR Status
# Rotterdam     == Int J Cancer 2013        == GSE41408    == HuEx 1.0 ST              == BCR/Metastasis Status
# Mayo Clinic I == Plos One 2013 (Decipher) == GSE46691    == HuEx 1.0 ST              == Metastasis Status
# GSE51066      == Prostate Cancer Prostatic Dis 2014 == GSE51066 == HuEx 1.0 ST       == Metastasis Status
# GSE37199      == Lancet Oncol 2012        == GSE37199    == HG-U133 Plus 2.0         == Prognosis Status (Good; Advanced Castration Resistant)

#######################

traits <- c('sample_id','patient_id','tissue','batch','platform','sample_type','age_at_diagnosis','ethnicity','race',
            'clinical_stage','clinical_t_stage','clinical_n_stage','clinical_m_stage',
            'pathological_stage','pathological_t_stage','pathological_n_stage','pathological_m_stage',
            'preop_psa','gleason_primary_pattern','gleason_secondary_pattern','gleason_tertiary_pattern',
            'gleason_group','gleason_score','time_to_death','os_status','time_to_bcr','bcr_status',
            'time_to_metastasis','metastasis_status','risk_group','treatment','additional_treatment')

### TCGA

exprData <- readRDS('data/rData/Expression_LogCPM_All_Genes_TCGA_PRAD.RDS')
dim(exprData)
exprData

meta.rna <- readRDS('data/rData/Metadata_RNAseq_TCGA_PRAD.RDS')
meta.rna
rownames(meta.rna) <- substr(rownames(meta.rna),start = 1, stop = 15)

pheno1 <- readRDS('data/rData/Clinical_TCGA_PRAD_11262019.RDS')
rownames(pheno1) <- paste0(rownames(pheno1), '-01')

pheno1 <- pheno1[match(rownames(meta.rna), rownames(pheno1)),]
rownames(pheno1) <- rownames(meta.rna)


pheno2 <- readRDS('data/rData/Clinical_TCGA_PRAD_With_PreopPSA_and_BCR.RDS')
pheno2 <- pheno2[match(rownames(meta.rna), rownames(pheno2)),]
pheno2

rownames(pheno2) <- rownames(meta.rna)

View(pheno1)
View(pheno2)


phenoData <- data.frame(matrix(NA, nrow=nrow(pheno1), ncol=length(traits)), stringsAsFactors = F)
phenoData

colnames(phenoData) <- traits
rownames(phenoData) <- rownames(pheno1)
colnames(phenoData)

phenoData$sample_id <- rownames(pheno1)
phenoData$patient_id <- substr(phenoData$sample_id, start = 1, stop = 12)
phenoData$platform <- 'Illumina'
phenoData$sample_type <- ifelse(grepl('-11', phenoData$sample_id), 'Normal', 'Primary')
phenoData$age_at_diagnosis <- pheno1$age_at_initial_pathologic_diagnosis
phenoData$ethnicity <- pheno1$ethnicity
phenoData$race <- pheno1$race
phenoData$clinical_t_stage <- pheno1$clinical_T
phenoData$clinical_n_stage <- pheno1$clinical_N
phenoData$clinical_m_stage <- pheno1$clinical_M
phenoData$pathological_t_stage <- pheno1$pathologic_T
phenoData$pathological_n_stage <- pheno1$pathologic_N
phenoData$pathological_m_stage <- pheno1$pathologic_M
phenoData$preop_psa <- pheno2$preop_psa
phenoData$gleason_primary_pattern <- pheno2$primary_pattern
phenoData$gleason_secondary_pattern <- pheno2$secondary_pattern
phenoData$gleason_tertiary_pattern <- pheno2$tertiary_pattern
phenoData$gleason_score <- phenoData$gleason_primary_pattern + phenoData$gleason_secondary_pattern
phenoData$gleason_group <- ifelse(is.na(phenoData$gleason_score), NA, paste(phenoData$gleason_primary_pattern, phenoData$gleason_secondary_pattern, sep='+'))

phenoData$time_to_death <- ifelse(!is.na(pheno2$days_to_death), 
                                  round(as.numeric(pheno2$days_to_death/365*12),2), 
                                  round(as.numeric(pheno2$days_to_last_followup/365*12),2))
phenoData$os_status <- ifelse(!is.na(pheno2$days_to_death), 1, 0)

phenoData$time_to_bcr <- ifelse(!is.na(pheno2$days_to_first_biochemical_recurrence), 
                                round(as.numeric(pheno2$days_to_first_biochemical_recurrence/365*12),2), 
                                round(as.numeric(pheno2$days_to_last_followup/365*12),2))
phenoData$bcr_status <- ifelse(!is.na(pheno2$days_to_first_biochemical_recurrence), 1, 0)


rownames(phenoData) == colnames(exprData)


eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/Database/TCGA_PRAD_eSet.RDS'))



#### CPC-GENE == Nature_2017 == GSE107299 == HuGene 2.0 ST & HTA 2.0 == Time to BCR (Download merged data from GEO)

eSet <- readRDS('data/rData/GSE107299_Merge_Prcessed_From_GEO_eSet.RDS')
exprData <- exprs(eSet)
keep <- which(startsWith(rownames(exprData),'ENSG'))
keep

exprData <- exprData[keep,]

pheno <- pData(eSet)
pheno

phenoData <- data.frame(matrix(NA, nrow=nrow(pheno), ncol=length(traits)), stringsAsFactors = F)
phenoData

colnames(phenoData) <- traits
rownames(phenoData) <- rownames(pheno)

phenoData$sample_id <- pheno$Sample.ID
phenoData$tissue <- pheno$tissue
phenoData$batch <- pheno$batch
phenoData$sample_type <- 'Primary'
phenoData$platform <- ifelse(phenoData$batch==2,'Affymetrix Human Transcriptome Array 2.0','Affymetrix Human Gene 2.0 ST Array')
phenoData$age_at_diagnosis <- pheno$Diagnosis.Age
phenoData$ethnicity <- pheno$Ethnicity.Category
phenoData$race <- pheno$Race.Category
phenoData$clinical_t_stage <- pheno$Clinical.T.Stage
phenoData$clinical_m_stage <- pheno$Clinical.M.Stage
phenoData$preop_psa <- pheno$Pre_Treatment_PSA
phenoData$gleason_primary_pattern <- pheno$Gleason.pattern.primary
phenoData$gleason_secondary_pattern <- pheno$Gleason.pattern.secondary
phenoData$gleason_tertiary_pattern <- pheno$Gleason.pattern.tertiary
phenoData$gleason_score <- phenoData$gleason_primary_pattern + phenoData$gleason_secondary_pattern
phenoData$gleason_group <- ifelse(is.na(phenoData$gleason_score), NA, paste(phenoData$gleason_primary_pattern, phenoData$gleason_secondary_pattern, sep='+'))
phenoData$time_to_bcr <- round(pheno$Time_to_Biochemical_Recurrence_or_Max_Follow_Up/365*12,2)
phenoData$bcr_status <- ifelse(pheno$Biochemical_Recurrence=='YES',1,0)
phenoData$treatment <- pheno$Treatment


rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/Database/GSE107299_eSet.RDS'))



### MSKCC ==_Cancer_Cell_2010 == GSE21034 == HuEx 1.0 ST == Time to BCR

eSet <- readRDS('data/rData/GSE21034_CDF24_GENCODE_RMA_eSet.RDS')
exprData <- exprs(eSet)
keep <- which(startsWith(rownames(exprData),'ENSG'))
keep

exprData <- exprData[keep,]
dim(exprData)

pheno <- pData(eSet)
pheno
View(pheno)

phenoData <- data.frame(matrix(NA, nrow=nrow(pheno), ncol=length(traits)), stringsAsFactors = F)
phenoData

colnames(phenoData) <- traits
rownames(phenoData) <- rownames(pheno)

phenoData$sample_id <- pheno$Sample.ID
phenoData$patient_id <- pheno$Patient.ID
phenoData$tissue <- 'Frozen'
phenoData$sample_type <- ifelse(pheno$Sample.Class=='Cell line', 'Cell line', pheno$Sample.Type)
phenoData$sample_type <- ifelse(is.na(phenoData$sample_type), 'Normal', phenoData$sample_type)
phenoData$platform <- 'Affymetrix Human Exon 1.0 ST Array'
phenoData$clinical_t_stage <- pheno$Neoplasm.American.Joint.Committee.on.Cancer.Clinical.Primary.Tumor.T.Stage
phenoData$gleason_primary_pattern <- pheno$Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer.1
phenoData$gleason_secondary_pattern <- pheno$Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer.2
phenoData$gleason_score <- phenoData$gleason_primary_pattern + phenoData$gleason_secondary_pattern
phenoData$gleason_group <- ifelse(is.na(phenoData$gleason_score), NA, paste(phenoData$gleason_primary_pattern, phenoData$gleason_secondary_pattern, sep='+'))
phenoData$time_to_bcr <- as.numeric(pheno$Disease.Free..Months.)
phenoData$bcr_status <- ifelse(pheno$Disease.Free.Status=='Recurred',1,0)

rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/Database/GSE21034_eSet.RDS'))



### DKFZ == Cancer Cell 2018 == cBioPortal == RNAseq: RPKM == Time to BCR  # one patient may be sequenced for multiple times

eSet <- readRDS('data/rData/DKFZ2018_cBioPortal_RPKM_eSet.RDS')
exprData <- exprs(eSet)
keep <- which(startsWith(rownames(exprData),'ENSG'))
keep

exprData <- exprData[keep,]
dim(exprData)

pheno <- pData(eSet)
pheno
View(pheno)

phenoData <- data.frame(matrix(NA, nrow=nrow(pheno), ncol=length(traits)), stringsAsFactors = F)
phenoData

colnames(phenoData) <- traits
rownames(phenoData) <- rownames(pheno)

colnames(phenoData)
phenoData$sample_id <- pheno$Sample.ID
phenoData$patient_id <- pheno$Patient.ID
phenoData$platform <- 'Illumina HiSeq 2000 (50bp paired-end)'
phenoData$sample_type <- 'Primary'
phenoData$age_at_diagnosis <- pheno$Diagnosis.Age
phenoData$pathological_stage <- pheno$Stage
phenoData$preop_psa <- pheno$Preop.PSA
phenoData$gleason_primary_pattern <- as.numeric(lapply(pheno$Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer, function(x) strsplit(x, '+', fixed=T)[[1]][1]))
phenoData$gleason_secondary_pattern <- as.numeric(lapply(pheno$Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer, function(x) strsplit(x, '+', fixed=T)[[1]][2]))
phenoData$gleason_score <- phenoData$gleason_primary_pattern + phenoData$gleason_secondary_pattern
phenoData$gleason_group <- ifelse(is.na(phenoData$gleason_score), NA, paste(phenoData$gleason_primary_pattern, phenoData$gleason_secondary_pattern, sep='+'))
phenoData$time_to_bcr <- as.numeric(pheno$Time.from.Surgery.to.BCR.Last.Follow.Up)
phenoData$bcr_status <- as.numeric(pheno$BCR.Status)

rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/Database/DKFZ2018_eSet.RDS'))



# GSE54460 == Cancer Research 2014 == GSE54460 == RNAseq: CPM == Time to BCR (Realignment and quantification)

eSet <- readRDS('data/rData/GSE54460_LogCPM_All_Genes_eSet.RDS')
exprData <- exprs(eSet)
keep <- which(startsWith(rownames(exprData),'ENSG'))
keep

exprData <- exprData[keep,]
dim(exprData)

pheno <- pData(eSet)
pheno
View(pheno)

phenoData <- data.frame(matrix(NA, nrow=nrow(pheno), ncol=length(traits)), stringsAsFactors = F)
phenoData

colnames(phenoData) <- traits
rownames(phenoData) <- rownames(pheno)
colnames(phenoData)

phenoData$sample_id <- pheno$title
phenoData$patient_id <- pheno$external_id
phenoData$tissue <- pheno$tissue
phenoData$platform <- 'Illumina HiSeq 2000 (50bp paired-end)'
phenoData$sample_type <- 'Primary'
phenoData$race <- pheno$race
phenoData$pathological_t_stage <- pheno$pathological_t_stage
phenoData$preop_psa <- as.numeric(pheno$preop_psa)
phenoData$gleason_primary_pattern <-as.numeric(substr(pheno$maybe_gleason_grade, 1, 1))
phenoData$gleason_secondary_pattern <- as.numeric(substr(pheno$maybe_gleason_grade, 2, 2))
phenoData$gleason_score <- phenoData$gleason_primary_pattern + phenoData$gleason_secondary_pattern
phenoData$gleason_group <- ifelse(is.na(phenoData$gleason_score), NA, paste(phenoData$gleason_primary_pattern, phenoData$gleason_secondary_pattern, sep='+'))
phenoData$time_to_bcr <- ifelse(pheno$bcr_status=='0', as.numeric(pheno$months_to_last_follow_up), as.numeric(pheno$months_to_bcr))
phenoData$bcr_status <- as.numeric(pheno$bcr_status)

phenoData$no_uniquely_mapped_reads <- pheno$no_of_uniquely_mapped_reads
phenoData$filter <- NA

samples <- c('SRR1164841','SRR1164843','SRR1164846','SRR1164847','SRR1164849','SRR1164851')
phenoData[samples,]$filter <- 'Duplicate'


rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/Database/Primary/GSE54460_eSet.RDS'))



# CamCap == Nature Genetics 2016 == GSE70768 == Illumina HumanHT-12 V4.0 == Time to BCR

eSet <- readRDS('data/rData/GSE70768_Expression_from_GEO_eSet.RDS')
exprData <- exprs(eSet)
keep <- which(startsWith(rownames(exprData),'ENSG'))
keep

exprData <- exprData[keep,]
dim(exprData)

pheno <- pData(eSet)
pheno
View(pheno)

phenoData <- data.frame(matrix(NA, nrow=nrow(pheno), ncol=length(traits)), stringsAsFactors = F)
phenoData

colnames(phenoData) <- traits
rownames(phenoData) <- rownames(pheno)
colnames(phenoData)

phenoData$sample_id <- pheno$description
phenoData$tissue <- pheno$tissue
phenoData$platform <- 'Illumina HumanHT-12 V4.0'
phenoData$sample_type <- ifelse(pheno$sample_type=='Tumour', 'Primary', 'Normal')
phenoData$age_at_diagnosis <- pheno$age
phenoData$pathological_stage <- pheno$pathological_stage
phenoData$clinical_stage <- pheno$clinical_stage
phenoData$preop_psa <- as.numeric(pheno$preop_psa)
phenoData$gleason_primary_pattern <-as.numeric(lapply(pheno$gleason_score, function(x) strsplit(x, '\\+|\\=')[[1]][2]))
phenoData$gleason_secondary_pattern <- as.numeric(lapply(pheno$gleason_score, function(x) strsplit(x, '\\+|\\=')[[1]][3]))
phenoData$gleason_score <- phenoData$gleason_primary_pattern + phenoData$gleason_secondary_pattern
phenoData$gleason_group <- ifelse(is.na(phenoData$gleason_score), NA, paste(phenoData$gleason_primary_pattern, phenoData$gleason_secondary_pattern, sep='+'))
phenoData$time_to_bcr <- ifelse(pheno$bcr_status=='Y', round(as.numeric(pheno$months_to_bcr),2), round(as.numeric(pheno$months_to_last_follow_up),2))
phenoData$bcr_status <- pheno$bcr_status
phenoData$bcr_status[phenoData$bcr_status=='N/A'] <- NA
phenoData$bcr_status <- ifelse(phenoData$bcr_status=='Y', 1, 0)


rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/Database/GSE70768_eSet.RDS'))




# Stockholm == Nature Genetics 2016 == GSE70769 == Illumina HumanHT-12 V4.0 == Time to BCR

eSet <- readRDS('data/rData/GSE70769_Expression_from_GEO_eSet.RDS')
exprData <- exprs(eSet)
keep <- which(startsWith(rownames(exprData),'ENSG'))
keep

exprData <- exprData[keep,]
dim(exprData)

pheno <- pData(eSet)
pheno
View(pheno)

phenoData <- data.frame(matrix(NA, nrow=nrow(pheno), ncol=length(traits)), stringsAsFactors = F)
phenoData

colnames(phenoData) <- traits
rownames(phenoData) <- rownames(pheno)
colnames(phenoData)

phenoData$sample_id <- gsub('tumour tissue_robotic radical prostatetctomy_','',pheno$title)
phenoData$tissue <- pheno$tissue
phenoData$platform <- 'Illumina HumanHT-12 V4.0'
phenoData$sample_type <- 'Primary'
phenoData$pathological_stage <- pheno$pathological_stage
phenoData$clinical_stage <- pheno$clinical_stage
phenoData$preop_psa <- as.numeric(pheno$preop_psa)
phenoData$gleason_primary_pattern <-as.numeric(lapply(pheno$gleason_score, function(x) strsplit(x, '\\+|\\=')[[1]][2]))
phenoData$gleason_secondary_pattern <- as.numeric(lapply(pheno$gleason_score, function(x) strsplit(x, '\\+|\\=')[[1]][3]))
phenoData$gleason_score <- phenoData$gleason_primary_pattern + phenoData$gleason_secondary_pattern
phenoData$gleason_group <- ifelse(is.na(phenoData$gleason_score), NA, paste(phenoData$gleason_primary_pattern, phenoData$gleason_secondary_pattern, sep='+'))
phenoData$time_to_bcr <- ifelse(pheno$bcr_status=='Y', round(as.numeric(pheno$months_to_bcr),2), round(as.numeric(pheno$months_to_last_follow_up),2))
phenoData$bcr_status <- pheno$bcr_status
phenoData$bcr_status[phenoData$bcr_status=='N/A'] <- NA
phenoData$bcr_status <- ifelse(phenoData$bcr_status=='Y', 1, 0)

rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/Database/GSE70769_eSet.RDS'))




# DESNT == Eur Urol Focus 2017 == GSE94767 == HuEx 1.0 ST == Time to BCR
eSet <- readRDS('data/rData/GSE94767_CDF24_GENCODE_RMA_eSet.RDS')
exprData <- exprs(eSet)
keep <- which(startsWith(rownames(exprData),'ENSG'))
keep

exprData <- exprData[keep,]
dim(exprData)

pheno <- pData(eSet)
pheno
View(pheno)

phenoData <- data.frame(matrix(NA, nrow=nrow(pheno), ncol=length(traits)), stringsAsFactors = F)
phenoData

colnames(phenoData) <- traits
rownames(phenoData) <- rownames(pheno)
colnames(phenoData)

phenoData$sample_id <- pheno$Sample.ID
phenoData$patient_id <- pheno$Donor.ID
phenoData$tissue <- 'Prostatectomy'
phenoData$batch <- pheno$Batch
phenoData$platform <- 'Affymetrix Human Exon 1.0 ST Array'
phenoData$sample_type <- pheno$material_type
phenoData$pathological_stage <- paste0(pheno$Pathology_Stage, pheno$Pathology_sub_stage)
phenoData$pathological_stage[phenoData$pathological_stage=='NANA'] <- NA
phenoData$preop_psa <- as.numeric(pheno$PSA_pre_prostatectomy)
phenoData$gleason_primary_pattern <-as.numeric(lapply(pheno$Gleason_Score, function(x) strsplit(x, '\\+|\\=')[[1]][1]))
phenoData$gleason_secondary_pattern <- as.numeric(lapply(pheno$Gleason_Score, function(x) strsplit(x, '\\+|\\=')[[1]][2]))
phenoData$gleason_score <- phenoData$gleason_primary_pattern + phenoData$gleason_secondary_pattern
phenoData$gleason_group <- ifelse(is.na(phenoData$gleason_score), NA, paste(phenoData$gleason_primary_pattern, phenoData$gleason_secondary_pattern, sep='+'))
phenoData$time_to_bcr <- as.numeric(pheno$BCR_FreeTime_months)
phenoData$bcr_status <- ifelse(pheno$BCR_Event=='TRUE', 1, 0)

rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/Database/GSE94767_eSet.RDS'))


# E-MTAB-6128   == Annuals of Oncology 2018 == E-MTAB-6128 == HuGene 2.0 ST            == Time to BCR

eSet <- readRDS('data/rData/E_MTAB_6128_CDF24_GENCODE_RMA_eSet.RDS')
exprData <- exprs(eSet)
keep <- which(startsWith(rownames(exprData),'ENSG'))
keep

exprData <- exprData[keep,]
dim(exprData)

pheno <- pData(eSet)
pheno
View(pheno)

phenoData <- data.frame(matrix(NA, nrow=nrow(pheno), ncol=length(traits)), stringsAsFactors = F)
phenoData

colnames(phenoData) <- traits
rownames(phenoData) <- rownames(pheno)
colnames(phenoData)

phenoData$sample_id <- pheno$Source.Name
phenoData$patient_id <- pheno$Characteristics.patient.id.
phenoData$platform <- 'Affymetrix Human Gene 2.0 ST Array'
phenoData$sample_type <- ifelse(grepl('normal',pheno$Characteristics.sampling.site.), 'Normal', 'Primary')
phenoData$age_at_diagnosis <- as.numeric(pheno$Characteristics.age.)
phenoData$preop_psa <- as.numeric(pheno$Characteristics.prostate.specific.antigen.measurement.)
phenoData$gleason_primary_pattern <-as.numeric(lapply(gsub('gleason score ','',pheno$Characteristics.gleason.score.), function(x) strsplit(x, '\\+|\\=')[[1]][1]))
phenoData$gleason_secondary_pattern <- as.numeric(lapply(gsub('gleason score ','',pheno$Characteristics.gleason.score.), function(x) strsplit(x, '\\+|\\=')[[1]][2]))
phenoData$gleason_score <- phenoData$gleason_primary_pattern + phenoData$gleason_secondary_pattern
phenoData$gleason_group <- ifelse(is.na(phenoData$gleason_score), NA, paste(phenoData$gleason_primary_pattern, phenoData$gleason_secondary_pattern, sep='+'))
phenoData$time_to_bcr <- as.numeric(pheno$Factor.Value.months.to.bcr.or.last.news.)
phenoData$bcr_status <- as.numeric(pheno$Factor.Value.bcr.event.)

rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/Database/E_MTAB_6128_eSet.RDS'))



# Belfast       == Journal of Clinical Oncology 2017  == GSE116918 == ADXPCv1a520642   == Time to BCR/Metastasis
eSet <- readRDS('data/rData/GSE116918_CDF24_GENCODE_RMA_eSet.RDS')
exprData <- exprs(eSet)
keep <- which(startsWith(rownames(exprData),'ENSG'))
keep

exprData <- exprData[keep,]
dim(exprData)

pheno <- pData(eSet)
pheno
View(pheno)

phenoData <- data.frame(matrix(NA, nrow=nrow(pheno), ncol=length(traits)), stringsAsFactors = F)
phenoData

colnames(phenoData) <- traits
rownames(phenoData) <- rownames(pheno)
colnames(phenoData)

phenoData$sample_id <- pheno$title
phenoData$platform <- '[ADXPCv1a520642] Almac Diagnostics Prostate Disease Specific Array (DSA)'
phenoData$sample_type <- 'Primary'
phenoData$age_at_diagnosis <- as.numeric(pheno$age)
phenoData$preop_psa <- as.numeric(pheno$preop_psa)
phenoData$gleason_score <- as.numeric(pheno$gleason_score)
phenoData$time_to_bcr <- as.numeric(pheno$months_to_bcr)
phenoData$bcr_status <- as.numeric(pheno$bcr_status)
phenoData$time_to_metastasis <- as.numeric(pheno$months_to_met)
phenoData$metastasis_status <- as.numeric(pheno$met_status)

rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/Database/GSE116918_eSet.RDS'))


# GSE25136      == Prostate 2019            == GSE25136    == HG-U133A                 == BCR Status

eSet <- readRDS('data/rData/GSE25136_CDF24_GENCODE_RMA_eSet.RDS')
exprData <- exprs(eSet)
keep <- which(startsWith(rownames(exprData),'ENSG'))
keep

exprData <- exprData[keep,]
dim(exprData)

pheno <- pData(eSet)
pheno
View(pheno)

phenoData <- data.frame(matrix(NA, nrow=nrow(pheno), ncol=length(traits)), stringsAsFactors = F)
phenoData

colnames(phenoData) <- traits
rownames(phenoData) <- rownames(pheno)
colnames(phenoData)

phenoData$sample_id <- pheno$description
phenoData$platform <- 'Affymetrix Human Genome U133A Array'
phenoData$sample_type <- 'Primary'
phenoData$bcr_status <- ifelse(pheno$recurrence_status=='Recurrent',1,0)

rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/Database/GSE25136_eSet.RDS'))


# GSE44353      == PNAS 2013 == GSE44353    == Prostate Cancer DASL Panel 1.5K         == BCR Status (data download from GEO)

eSet <- readRDS('data/rData/GSE44353_from_GEO_eSet.RDS')
exprData <- exprs(eSet)
keep <- which(startsWith(rownames(exprData),'ENSG'))
keep

exprData <- exprData[keep,]
dim(exprData)

pheno <- pData(eSet)
pheno
View(pheno)

phenoData <- data.frame(matrix(NA, nrow=nrow(pheno), ncol=length(traits)), stringsAsFactors = F)
phenoData

colnames(phenoData) <- traits
rownames(phenoData) <- rownames(pheno)
colnames(phenoData)

phenoData$sample_id <- pheno$title
phenoData$platform <- 'Illumina Custom Huamn Prostate Cancer DASL Panel 1.5K'
phenoData$sample_type <- 'Primary'
phenoData$pathological_t_stage <- pheno$pathological_t_stage
phenoData$gleason_score <- as.numeric(pheno$gleason_score)
phenoData$bcr_status <- as.numeric(pheno$bcr_status)

rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/Database/GSE44353_eSet.RDS'))


# Rotterdam     == Int J Cancer 2013        == GSE41408    == HuEx 1.0 ST              == BCR/Metastasis Status

eSet <- readRDS('data/rData/GSE41408_CDF24_GENCODE_RMA_eSet.RDS')
exprData <- exprs(eSet)
keep <- which(startsWith(rownames(exprData),'ENSG'))
keep

exprData <- exprData[keep,]
dim(exprData)

pheno <- pData(eSet)
pheno
View(pheno)

phenoData <- data.frame(matrix(NA, nrow=nrow(pheno), ncol=length(traits)), stringsAsFactors = F)
phenoData

colnames(phenoData) <- traits
rownames(phenoData) <- rownames(pheno)
colnames(phenoData)

phenoData$sample_id <- pheno$title
phenoData$tissue <- pheno$tissue
phenoData$platform <- 'Affymetrix Human Exon 1.0 ST Array'
phenoData$sample_type <- 'Primary'
phenoData$pathological_t_stage <- pheno$pathological_t_stage
phenoData$preop_psa <- as.numeric(phenoData$preop_psa)
phenoData$gleason_primary_pattern <- as.numeric(lapply(pheno$gleason_score, function(x) strsplit(x, '\\+|\\=')[[1]][1]))
phenoData$gleason_secondary_pattern <- as.numeric(lapply(pheno$gleason_score, function(x) strsplit(x, '\\+|\\=')[[1]][2]))
phenoData$gleason_score <- phenoData$gleason_primary_pattern + phenoData$gleason_secondary_pattern
phenoData$gleason_group <- ifelse(is.na(phenoData$gleason_score), NA, paste(phenoData$gleason_primary_pattern, phenoData$gleason_secondary_pattern, sep='+'))

phenoData$bcr_status <- ifelse(pheno$bcr_status=='yes', 1, 0)
phenoData$metastasis_status <- ifelse(pheno$met_status=='yes', 1, 0)

rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/Database/GSE41408_eSet.RDS'))



# Mayo Clinic I == Plos One 2013 (Decipher) == GSE46691    == HuEx 1.0 ST              == Metastasis Status

eSet <- readRDS('data/rData/GSE46691_CDF24_GENCODE_RMA_eSet.RDS')
exprData <- exprs(eSet)
keep <- which(startsWith(rownames(exprData),'ENSG'))
keep

exprData <- exprData[keep,]
dim(exprData)

pheno <- pData(eSet)
pheno
View(pheno)


phenoData <- data.frame(matrix(NA, nrow=nrow(pheno), ncol=length(traits)), stringsAsFactors = F)
phenoData

colnames(phenoData) <- traits
rownames(phenoData) <- rownames(pheno)
colnames(phenoData)

phenoData$sample_id <- pheno$description
phenoData$platform <- 'Affymetrix Human Exon 1.0 ST Array'
phenoData$sample_type <- 'Primary'
phenoData$gleason_score <- as.numeric(pheno$gleason_score)
phenoData$metastasis_status <- as.numeric(pheno$metastatic_event)

rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/Database/GSE46691_eSet.RDS'))


# GSE51066      == Prostate Cancer Prostatic Dis 2014 == GSE51066 == HuEx 1.0 ST       == Metastasis Status

eSet <- readRDS('data/rData/GSE51066_CDF24_GENCODE_RMA_eSet.RDS')
exprData <- exprs(eSet)
keep <- which(startsWith(rownames(exprData),'ENSG'))
keep

exprData <- exprData[keep,]
dim(exprData)

pheno <- pData(eSet)
pheno
View(pheno)


phenoData <- data.frame(matrix(NA, nrow=nrow(pheno), ncol=length(traits)), stringsAsFactors = F)
phenoData

colnames(phenoData) <- traits
rownames(phenoData) <- rownames(pheno)
colnames(phenoData)

phenoData$sample_id <- pheno$title
phenoData$tissue <- 'FFPE'
phenoData$platform <- 'Affymetrix Human Exon 1.0 ST Array'
phenoData$sample_type <- 'Primary'
phenoData$metastasis_status <- as.numeric(pheno$metastatic_event)

rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/Database/GSE51066_eSet.RDS'))


# GSE37199      == Lancet Oncol 2012        == GSE37199    == HG-U133 Plus 2.0         == Prognosis Status (Good; Advanced Castration Resistant)

eSet <- readRDS('data/rData/GSE37199_CDF24_GENCODE_RMA_eSet.RDS')
exprData <- exprs(eSet)
keep <- which(startsWith(rownames(exprData),'ENSG'))
keep

exprData <- exprData[keep,]
dim(exprData)

pheno <- pData(eSet)
pheno
View(pheno)


phenoData <- data.frame(matrix(NA, nrow=nrow(pheno), ncol=length(traits)), stringsAsFactors = F)
phenoData

colnames(phenoData) <- traits
rownames(phenoData) <- rownames(pheno)
colnames(phenoData)

phenoData$sample_id <- pheno$title
phenoData$patient_id <- pheno$patient
phenoData$tissue <- pheno$tissue
phenoData$platform <- 'Affymetrix Human Genome U133 Plus 2.0 Array'
phenoData$sample_type <- 'Primary'
phenoData$risk_group <- pheno$prognosis_status

rownames(phenoData) == colnames(exprData)

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file=paste0('data/Database/GSE37199_eSet.RDS'))


########################################################################################



