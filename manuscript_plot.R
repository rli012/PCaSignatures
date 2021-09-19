setwd('~/bigdata/PCa/')

library(readxl)
library(stringi)
library(VennDiagram)
library(Vennerable)
library(corrplot)

### KEGG Enrichment of 1032 genes

kegg <- read.csv('data/KEGG_All_Signature_Genes.csv', header = T, stringsAsFactors = F)

EnrichmentBubblePlotFun <- function(dataForBubblePlot) {
  
  p <- ggplot(dataForBubblePlot, mapping=aes(x=Description, y=Fold.Enrichment, #y=-log10(Benjamini), #y=Fold.Enrichment
                                             color=BH.Adj.P,size=Count)) +
    geom_point()+ coord_flip() +
    scale_x_discrete(limits=rev(unique(dataForBubblePlot$Description))) +
    #scale_x_discrete(limits=Order)+
    scale_colour_gradientn(limits=c(0,0.05),
                           colors= c("red","yellow","green")) + #
    #facet_wrap(~Comparison) +
    #facet_grid(Regulation~Comparison) + # scales=free
    xlab('')+ylab('Fold Enrichment') + #ggtitle("") + 
    guides(shape = guide_legend(order=1),
           colour = guide_colourbar(order=2, title = 'FDR')) + #'P Value\n(Benjamini)'))
    theme_bw()+theme(axis.line = element_line(colour = "black"),
                     panel.grid.minor = element_blank(),
                     panel.border = element_rect(colour='black'),
                     panel.background = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5, size=20)) +
    theme(axis.text=element_text(size=14, color='black', face = 'bold'),
          axis.text.x =element_text(size=14, color='black', face = 'bold', angle=0, hjust=0.5),
          axis.title=element_text(size=16, face = 'bold')) +
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    theme(#strip.text = element_text(size = 14),
      legend.key.size = unit(0.8,'cm'))
  
  
  return (p)
  
}

r1 <- sapply(kegg$Count.List.Total, function(x) as.numeric(strsplit(x, '/')[[1]][1])/as.numeric(strsplit(x, '/')[[1]][2]))
r2 <- sapply(kegg$Pop.Hits.Pop.Total, function(x) as.numeric(strsplit(x, '/')[[1]][1])/as.numeric(strsplit(x, '/')[[1]][2]))
kegg$Fold.Enrichment <- r1/r2
kegg$Count <- sapply(kegg$Count.List.Total, function(x) as.numeric(strsplit(x, '/')[[1]][1]))
p <- EnrichmentBubblePlotFun(kegg)
p

topptx(p, filename = 'figures/Figure_KEGG_1032_Genes.pptx', width=10, height = 7)

### signatures genes in pathways
signatures <- readxl::read_excel('../PCaDB/www/downloads/PCaDB_Prognostic_Signatures_Gene_List.xlsx')

signature.names <- unique(signatures$Signature)
signature.names

mat <- data.frame(matrix(rep(0, nrow(kegg)*length(signature.names)), 
                         nrow=nrow(kegg), ncol=length(signature.names)),
                  stringsAsFactors = F)

rownames(mat) <- kegg$Description
colnames(mat) <- signature.names

genes <- c()
for (i in 1:nrow(kegg)) {
  pathway <- kegg$Description[i]
  symbols <- strsplit(kegg$Symbol[i], '; ')[[1]]
  genes <- c(genes, symbols)
  
  for (sig in signature.names) {
    idx <- which(signatures$Signature==sig)
    sig.genes <- signatures$HGNC.Symbol[idx]
    
    ovlp <- intersect(sig.genes, symbols)
    
    mat[pathway, sig] <- length(ovlp)
  }
  
  
}

View(mat)

c(mat)


dataForBubblePlot <- data.frame(Count=unlist(c(mat)), 
                                Signature=rep(signature.names,each=nrow(kegg)),
                                Pathway=rep(kegg$Description, length(signature.names)),
                                stringsAsFactors = F)

saveRDS(dataForBubblePlot, file='manuscript/Figure_KEGG_Signature_Overlap_dataForBubblePlot.RDS')


#View(dataForBubblePlot)

dataForBubblePlot$Signature <- factor(dataForBubblePlot$Signature, levels = signature.names)
dataForBubblePlot$Pathway <- factor(dataForBubblePlot$Pathway, levels = rev(kegg$Description))

dataForBubblePlot$Count[dataForBubblePlot$Count==0] <- NA

p <- ggplot(dataForBubblePlot, mapping=aes(x=as.numeric(Signature)+0.5, y=as.numeric(Pathway)+0.5)) +
  geom_point(aes(color='red', size=Count), color= google.red, alpha=1)+ 
  scale_x_continuous(minor_breaks = seq(1, 30, 1), limits = c(1,31),expand = c(0, 0), position = "top",
                     breaks = seq(1,30,1), labels = unique(dataForBubblePlot$Signature)) +
  scale_y_continuous(minor_breaks = seq(1, 34, 1), limits = c(1,35),expand = c(0, 0),
                     breaks = seq(1,34,1), labels = rev(unique(dataForBubblePlot$Pathway))) +
  # geom_text(data=dataForBubbleAnnot, mapping=aes(x=as.numeric(Signature)+0.5, 
  #                                                y=as.numeric(Pathway)+0.5, 
  #                                                label=Count), size=4.5, color=google.blue) +
  theme_bw() +
  
  theme(panel.grid.minor = element_line(colour="gray", size=0.5),
        panel.grid.major = element_line(colour="gray", size=0.5),
        panel.border = element_rect(colour = "gray", fill=NA, size = 1),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust=-1.5, hjust=0, size=14, face = 'bold'),
        axis.text.y = element_text(vjust = -0.5, size=14, face = 'bold'),
        axis.ticks = element_blank()) +
  theme(legend.position = 'bottom')

p

topptx(p, filename = 'Figure_KEGG_Signature_Overlaps.pptx', width=12, height = 9)


###
signatures <- readxl::read_excel('figures/Table_PCa_Prognosis_Signatures.xlsx',
                                 sheet = 'Signatures')
signatures

sort(table(signatures$`HGNC Symbol`), decreasing = T)

signature.name <- unique(signatures$Signature)

signature.list <- list()

for (sig in signature.name) {
  
  genes <- signatures$`Ensembl ID`[signatures$Signature==sig]
  
  if (sum(is.na(genes))>0) {
    genes[is.na(genes)] <- stringi::stri_rand_strings(2, 16)
  }
  
  signature.list[[sig]] <- genes
  
  
}


ovlp <- outer(signature.list, signature.list, Vectorize(function(x, y) sum(x %in% y)))
View(ovlp)


dataForBubblePlot <- data.frame(Count=c(ovlp), Signature1=rep(paste0(unique(signatures$Signature),' (', unlist(lapply(signature.list, length)), ')'),each=30),
                                Signature2=rep(paste0(unique(signatures$Signature),' (', unlist(lapply(signature.list, length)), ')'),30),
                                stringsAsFactors = F)

saveRDS(dataForBubblePlot, file='manuscript/Figure_Signature_Overlap_dataForBubblePlot.RDS')

dataForBubblePlot <- readRDS(file='figures/Figure_Signature_Overlap_dataForBubblePlot.RDS')

#View(dataForBubblePlot)

dataForBubblePlot$Signature1 <- factor(dataForBubblePlot$Signature1, levels = unique(dataForBubblePlot$Signature1))
dataForBubblePlot$Signature2 <- factor(dataForBubblePlot$Signature2, levels = rev(unique(dataForBubblePlot$Signature2)))


dataForBubblePlot$Count[dataForBubblePlot$Count==0] <- NA
dataForBubbleAnnot <- dataForBubblePlot

idx <- which(as.character(dataForBubblePlot$Signature1) < as.character(dataForBubblePlot$Signature2))
dataForBubblePlot$Count[idx] <- NA


idx <- which(as.character(dataForBubbleAnnot$Signature1) >= as.character(dataForBubbleAnnot$Signature2))
dataForBubbleAnnot$Count[idx] <- NA



p <- ggplot(dataForBubblePlot, mapping=aes(x=as.numeric(Signature1)+0.5, y=as.numeric(Signature2)+0.5)) +
  geom_point(aes(color='red', size=Count), color= google.red, alpha=1)+ 
  scale_y_continuous(minor_breaks = seq(1, 30, 1), limits = c(1,31),expand = c(0, 0),
                     breaks = seq(1,30,1), labels = rev(unique(dataForBubblePlot$Signature1))) +
  scale_x_continuous(minor_breaks = seq(1, 30, 1), limits = c(1,31),expand = c(0, 0), position = "top",
                     breaks = seq(1,30,1), labels = unique(dataForBubblePlot$Signature2)) +
  geom_text(data=dataForBubbleAnnot, mapping=aes(x=as.numeric(Signature1)+0.5, 
                                                 y=as.numeric(Signature2)+0.5, 
                                                 label=Count), size=4.5, color=google.blue) +
  theme_bw() +
  
  theme(panel.grid.minor = element_line(colour="gray", size=0.5),
        panel.grid.major = element_line(colour="gray", size=0.5),
        panel.border = element_rect(colour = "gray", fill=NA, size = 1),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust=-1.5, hjust=0, size=14, face = 'bold'),
        axis.text.y = element_text(vjust = -0.5, size=14, face = 'bold'),
        axis.ticks = element_blank()) +
  theme(legend.position = 'none')

p
#library(eoffice)
topptx(p, filename = 'figures/Figure_Signature_Overlaps.pptx', width=8, height = 8)

#corrplot(ovlp, is.corr = FALSE, method = "circle", cl.lim = c(0, 250))

#################

## Circos

bed1 <- c()
bed2 <- c()

for (i in 1:29) {
  
  .g1 <- signature.list[[i]]
  
  for (j in (i+1):30) {
    
    .g2 <- signature.list[[j]]
    
    .ovlp <- intersect(.g1, .g2)
    
    if (length(.ovlp) > 0) {
          idx1 <- match(.ovlp, .g1)
          idx2 <- match(.ovlp, .g2)
          
          bed1 <- rbind(bed1, cbind(rep(names(signature.list)[i], length(idx1)), idx1, idx1))
          bed2 <- rbind(bed2, cbind(rep(names(signature.list)[j], length(idx2)), idx2, idx2))
    }
    
  }
  
}



bed1
bed2

bed1 <- data.frame(bed1, stringsAsFactors = F)
colnames(bed1) <- c('chr','start','end')

bed1$start <- bed1$end <- as.numeric(bed1$start)

bed2 <- data.frame(bed2, stringsAsFactors = F)
colnames(bed2) <- c('chr','start','end')

bed2$start <- bed2$end <- as.numeric(bed2$start)



set.seed(999)
set.seed(777)
#n = nrow(bed1)
n <- 1237
df = data.frame(signatures = rep(names(unlist(lapply(signature.list, length))), unlist(lapply(signature.list, length))),
                x = as.numeric(unlist(sapply(unlist(lapply(signature.list, length)), function(v) 1:v))), 
                y = runif(n))






library("RColorBrewer")
display.brewer.all()
cols <- sample(brewer.pal(n = 10, name = "Set3"),10)
cols 
#  [1] "#80B1D3" "#8DD3C7" "#BEBADA" "#B3DE69" "#BC80BD" "#FCCDE5" "#D9D9D9" "#FFFFB3" "#FDB462"
# [10] "#FB8072"
bg.col <- rep(cols,3)
#bg.col <- rep(brewer.pal(n = 10, name = "Paired"),10)



library(circlize)

circos.par(cell.padding = c(0.02, 0.02, 0.02, 0.02), "track.height" = 0.1, circle.margin=c(0.25))# 
circos.initialize(df$signatures, x = df$x)

all_facing_options = unlist(lapply(signature.list, length))

circos.track(df$signatures, y = df$y, bg.col = bg.col,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter,
                           CELL_META$cell.ylim[2] + mm_y(1),
                           CELL_META$sector.index,
                           facing = 'clockwise', 
                           cex = 1,
                           adj = c(0,0.5))
               circos.text(CELL_META$xcenter, CELL_META$ycenter, 
                           all_facing_options[CELL_META$sector.numeric.index],
                           facing = "inside", niceFacing = TRUE, cex=1) #1
             #   circos.axis(labels.cex = 0.001)
             }
             )


circos.genomicLink(bed1, bed2, border = NA, col = brewer.pal(n = 12, name = "Set3")[6])

circos.clear()


###############################################################


gene.stats <- data.frame(sort(table(signatures$`Ensembl ID`), decreasing = T),
                         stringsAsFactors = F)

#which(is.na(gene.stats$`Ensembl ID`))

colnames(gene.stats) <- c('Ensembl ID', 'Count')


match(gene.stats$`Ensembl ID`[1], signatures$`Ensembl ID`)

.sig.count <- c()

for (.g in signatures$`Ensembl ID`) {
  .count <- gene.stats$Count[which(gene.stats$`Ensembl ID`==.g)]
  .sig <- paste(signatures$Signature[which(signatures$`Ensembl ID`==.g)], collapse  = ';')
  
  .sig.count <- rbind(.sig.count, c(.count, .sig))
  
}

.sig.count

signatures$Count <- as.numeric(.sig.count[,1])

signatures$Common <- as.character(.sig.count[,2])

View(signatures)

# saveRDS(signatures, file = 'manuscript/Signatures.RDS')



##########

gene.stats <- gene.stats[-which(gene.stats$`Ensembl ID`=='NA'),]
gene.stats

gene.stats$`HGNC Symbol` <- signatures$`HGNC Symbol`[match(gene.stats$`Ensembl ID`, signatures$`Ensembl ID`)]

write.table(gene.stats, file = 'manuscript/Gene_Stats.txt',
            sep = '\t', quote = F, row.names = F)



dataForHeatmap <- matrix(rep(0, nrow(gene.stats)*length(signature.list)),
                         nrow = nrow(gene.stats), ncol = length(signature.list))


dataForHeatmap <- data.frame(dataForHeatmap, stringsAsFactors = F)
colnames(dataForHeatmap) <- names(signature.list)
rownames(dataForHeatmap) <- signatures$`HGNC Symbol`[match(gene.stats$`Ensembl ID`, signatures$`Ensembl ID`)]


.sig.gene <- paste0(signatures$Signature, ': ', signatures$`HGNC Symbol`)

for (.s in colnames(dataForHeatmap)) {
  for (.g in rownames(dataForHeatmap)) {
    
    if (paste(.s, .g, sep=': ') %in% .sig.gene) {
      
      dataForHeatmap[.g, .s] <- 1
      
    }
  }
}

rownames(dataForHeatmap) <- paste0(signatures$`HGNC Symbol`[match(gene.stats$`Ensembl ID`, signatures$`Ensembl ID`)],
                                   ' (', gene.stats$Count, ')')

View(dataForHeatmap)

# saveRDS(dataForHeatmap, file = 'manuscript/Figure_Signature_Overlap_dataForHeatmap.RDS')

dataForHeatmap <- readRDS(file = 'figures/Figure_Signature_Overlap_dataForHeatmap.RDS')


dataForHeatmap <- dataForHeatmap[which(apply(dataForHeatmap, 1, sum)>=2),]
dataForHeatmap <- dataForHeatmap[which(apply(dataForHeatmap, 1, sum)>=3),]
dim(dataForHeatmap)

View(dataForHeatmap)

library(RColorBrewer)
library(ComplexHeatmap)

cols = brewer.pal(4, "Reds")
cols = colorRampPalette(cols)(10)
cols
# col_fun = colorRampPalette(rev(c(cols[10],'white')), space = "Lab")(11)
# col_fun

library(circlize)
col_fun = colorRampPalette(rev(c(cols[10],cols[1])), space = "Lab")(2)

#col_fun = colorRampPalette(rev(c(cols[10],cols[1])), space = "Lab")(100)

rownames(dataForHeatmap) <- gsub('\\s\\(\\S+\\)', '', rownames(dataForHeatmap))

ht <- Heatmap(as.matrix(dataForHeatmap),
              #name = 'Expression',
              
              # COLOR
              #col = colorRampPalette(rev(c("red",'white','blue')), space = "Lab")(100),
              col=col_fun,
              na_col = 'gray',
              rect_gp = gpar(col = "grey", lwd = 0.5),
              
              # MAIN PANEL
              column_title = NULL,
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              show_column_dend = FALSE,
              show_row_names = TRUE,
              row_names_side = "left",
              column_names_side = 'top',
              show_column_names = TRUE,
              column_names_rot = 90,
              column_names_gp = gpar(fontsize = 14, fontface='bold'),
              row_names_gp = gpar(fontsize = 12, fontface='bold'),
              #column_names_max_height = unit(3, 'cm'),
              #column_split = factor(phenoData$Day,
              #                      levels=str_sort(unique(phenoData$Day), numeric = T)),
              
              #column_order = rownames(phenoData),
              
              # ANNOTATION
              #top_annotation = topAnnotation,
              
              # ADD TEXT
              # cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
              #   grid.text(annoMatrix[i, j], x, y, gp = gpar(fontsize = 12, col = "black", fill = 'black'))
              # },
              
              # LEGEND
              # heatmap_legend_param = list(
              #   at = c(-2,0,2,4,6,8,10,12),
              #   #labels = c("Negative",'Positive'),
              #   title = "",
              #   title_position = 'leftcenter-rot',
              #   #legend_height = unit(3, "cm"),
              #   adjust = c("right", "top")
              # ),
              show_heatmap_legend = FALSE
)

p <- draw(ht,annotation_legend_side = "right",row_dend_side = "left")

p

topptx(p, filename = 'figures/Figure_Heatmap_Common_Gene.pptx', width=7, height = 8)


#############################################################

###### Differential expression in TCGA

library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)

phenoData <- readRDS('~/bigdata/LABDATA/Collaboration/miR7/data/Metadata_RNAseq_TCGA_PRAD.RDS')

countsMatrix <- readRDS('~/bigdata/LABDATA/Collaboration/miR7/data/RNAseq_Counts_TCGA_PRAD.RDS')

annoData <- readRDS('~/bigdata/LABDATA/Collaboration/miR7/data/Homo_Sapiens_Gene_Annotation_ENSEMBL_HGNC_ENTREZ.RDS')


dge <-  DGEList(counts = countsMatrix)

### TMM normalization
dge = calcNormFactors(dge, method = 'TMM')

### Filter out low-expression genes (cpm>1 in at least 50% of the samples)
keep <- rowSums(cpm(dge) > 1) >= 0.5*ncol(countsMatrix)
sum(keep)
dge <- dge[keep,,keep.lib.sizes = TRUE]

### Voom normalization
v <- voom(dge, design=NULL, plot = FALSE)

exprAfterVoom <- v$E ### for visualization
exprLogCPM <- cpm(dge,log = TRUE) ### for visualization
exprLogCPM

### Prepare comparison matrix
group <- factor(phenoData$sample_type)
#group <- factor(phenoData$HistologyGrade)

design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design
contrast.matrix <- makeContrasts(contrasts='PrimaryTumor - SolidTissueNormal',
                                 levels=design)
contrast.matrix

### Differential gene expression analysis (limma)

fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


### Report DEGs
dgeTable <- topTable(fit2, coef=1, n=Inf, adjust.method='BH', sort.by='p')
#dgeTable <- topTable(fit2, coef='Moderate', n=Inf, adjust.method='BH', sort.by='p')
View(dgeTable)


### Map Ensembl ID to gene symbol
idx <- match(rownames(dgeTable), annoData$ensembl_id)
dgeTable$Symbol <- annoData$gene_name[idx]
dgeTable$Biotype <- annoData$gene_biotype[idx]


View(dgeTable)

dgeTable <- read.csv('figures/PCaDB_Differential_Expression_TCGA.csv', header = T, stringsAsFactors = F)
dgeTable

dataForVolcanoPlot <- dgeTable

dataForVolcanoPlot$Ensembl <- rownames(dataForVolcanoPlot)

logFcThreshold <- log2(2)
adjPvalThreshold <- 0.01

dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, 
                                     logFC < logFcThreshold | adj.P.Val > adjPvalThreshold)] <- 'NS'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, 
                                     logFC >= logFcThreshold & adj.P.Val <= adjPvalThreshold)] <- 'UP'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, 
                                     logFC <= -logFcThreshold & adj.P.Val <= adjPvalThreshold)] <- 'DOWN'

dataForVolcanoPlot$Significance.Label <- dataForVolcanoPlot$Significance

dataForVolcanoPlot$Significance <- 'NS'

signatures <- readRDS(file = 'figures/Signatures.RDS')

dataForVolcanoPlot$Signatures <- signatures$Common[match(dataForVolcanoPlot$Ensembl, signatures$`Ensembl ID`)]
dataForVolcanoPlot$Signature.Count <- signatures$Count[match(dataForVolcanoPlot$Ensembl, signatures$`Ensembl ID`)]


#saveRDS(dataForVolcanoPlot, file='manuscript/Figure_dataForVolcanoPlot_TCGA.RDS')


cols <- c('UP'=google.red, 'NS'='gray','DOWN'=google.green)
cols <- c('UP'=google.red, 'NS'='gray','DOWN'=google.blue)


dataForVolcanoPlot <- dataForVolcanoPlot[which(dataForVolcanoPlot$Signature.Count>=1),]

dim(dataForVolcanoPlot)
View(dataForVolcanoPlot)

p <- ggplot(dataForVolcanoPlot, aes(x = logFC, y = -log10(adj.P.Val))) +
  #xlim(-2,2) +
  labs(x=expression(bold('Log'['2']*'(Fold Change)')), 
       y=(expression(bold('-Log'['10']*'(FDR)'))), 
       title=NULL) +
  geom_point(color='gray', alpha=1, size=2) +
  geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold), 
             color='darkgreen',linetype='dashed')+
  xlim(-4,4) +
  #scale_x_continuous(breaks=c(-4,-2,0,2,4,6,8,10)) +
  #scale_y_continuous(expand = c(0.3, 0)) +
  #scale_color_manual(values = c('#4285F4',"gray", '#FBBC05')) +
  scale_color_manual(values = cols) +
  #facet_wrap(~Comparison, ncol = 2) +
  # geom_point(data = subset(dataForVolcanoPlot,
  #                          miRNA.ID %in% c('hsa-miR-7-1-3p','hsa-miR-7-5p')),
  #            aes(x=logFC, y=-log10(adj.P.Val)), size=2, color='blue', fill='gray',shape=21)+

  
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance.Label == 'UP' & Signature.Count>=1),
             aes(x=logFC, y=-log10(adj.P.Val)), fill=google.red, 
             size=2, color=google.red,shape=21)+
  
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance.Label == 'DOWN' & Signature.Count>=1),
             aes(x=logFC, y=-log10(adj.P.Val)), fill=google.blue, 
             size=2, color=google.blue,shape=21)+
  scale_fill_manual(values = cols) +
  
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance.Label == 'UP' & Signature.Count>=3),
                  aes(label = Symbol),nudge_x = 0.1,
                  segment.alpha = 0.5, #segment.size = 0.5,
                  min.segment.length = 0,
                  size = 4, color='black', segment.color = 'black') +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance.Label == 'DOWN' & Signature.Count>=3),
                  aes(label = Symbol),
                  segment.alpha = 0.5, #segment.size = 0.5,
                  #min.segment.length = 5,
                  size = 4, color='black', segment.color = 'black') +
theme_bw() +
  theme(axis.line = element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        strip.text = element_text(size=14, face='bold')) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


p

#write.table(dataForVolcanoPlot, file='manuscript/RNAseq_DEG_Tumor_vs_Normal.txt',
#            sep = '\t', quote = F, row.names = F)

topptx(p, filename = 'figures/Figure_Volcano_DEG_blue_red.pptx', width=6, height = 6)




#########################

eSet <- readRDS('data/Database/Primary/TCGA_PRAD_eSet.RDS')

phenoData <- pData(eSet)
exprData <- exprs(eSet)

dim(phenoData)

idx <- which(phenoData$sample_type=='Primary')
idx


phenoData <- phenoData[idx,]
exprData <- exprData[,idx]


library(survival)
library(survminer)

################## Survival Analysis


### RFS
time.to.event <- as.numeric(phenoData$time_to_bcr)
time.to.event

event.status <- phenoData$bcr_status
event.status

gene.stats$`Ensembl ID` <- as.character(gene.stats$`Ensembl ID`)
gene.stats$`HGNC Symbol` <- as.character(gene.stats$`HGNC Symbol`)

coxTable <- c()

for (.g in gene.stats$`Ensembl ID`) {
  
  if (!.g %in% rownames(exprData)) {
    coeffs <- rep(NA,5)
    coxTable <- rbind(coxTable, coeffs)
    next
  }
  
  expr <- exprData[.g,]
  coxtest <- coxph(Surv(time.to.event, event.status) ~ expr)
  summcph <- summary(coxtest)
  
  coeffs <- c(summcph$coefficients[,1:2], summcph$conf.int[,3:4], 
              summcph$coefficients[,5])
  
  
  coxTable <- rbind(coxTable, coeffs)
  
}

colnames(coxTable) <- c('Coef','HR','Lower95','Upper95','P')
rownames(coxTable) <- gene.stats$`Ensembl ID`

coxTable

coxTable <- data.frame(coxTable, stringsAsFactors = F)

coxTable$FDR <- p.adjust(coxTable$P, method = 'BH')

coxTable$Symbol <- gene.stats$`HGNC Symbol`

saveRDS(coxTable, file='manuscript/CoxPH_Genes_TCGA.RDS')

View(coxTable)


kmTable <- c()

for (.g in gene.stats$`Ensembl ID`) {
  
  if (!.g %in% rownames(exprData)) {
    coeffs <- rep(NA,5)
    kmTable <- rbind(kmTable, coeffs)
    next
  }
  
  expr <- exprData[.g,]
  risk.group <- expr > median(expr, na.rm = T)
  
  n.high <- sum(risk.group, na.rm=T)
  n.low <- sum(!risk.group, na.rm=T)
  
  sdf <- survdiff(Surv(time.to.event, event.status) ~ risk.group)
  p.val <- pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE)
  #p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  
  hr = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
  upper95 = exp(log(hr) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
  lower95 = exp(log(hr) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
  
  coeffs <- c(hr, lower95, upper95, p.val)
  
  kmTable <- rbind(kmTable, coeffs)
  
}

colnames(kmTable) <- c('HR','Lower95','Upper95','P')
rownames(kmTable) <- gene.stats$`Ensembl ID`

kmTable <- data.frame(kmTable, stringsAsFactors = F)

kmTable$FDR <- p.adjust(kmTable$P, method = 'BH')

kmTable$Symbol <- gene.stats$`HGNC Symbol`

saveRDS(kmTable, file='manuscript/KM_Genes_TCGA.RDS')


sum(coxTable$P[1:42] < 0.01, na.rm = T)

#############

kmTable <- readRDS(file='figures/CoxPH_1032_Genes_TCGA.RDS')

#dataForVolcanoPlot <- coxTable
dataForVolcanoPlot <- kmTable
dataForVolcanoPlot$Ensembl <- rownames(dataForVolcanoPlot)

hrThreshold <- log2(1)
adjPvalThreshold <- 0.05

dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, FDR > adjPvalThreshold)] <- 'NS'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, 
                                     log2(HR) >= hrThreshold & FDR <= adjPvalThreshold)] <- 'UP'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, 
                                     log2(HR) <= hrThreshold & FDR <= adjPvalThreshold)] <- 'DOWN'

dataForVolcanoPlot$Significance.Label <- dataForVolcanoPlot$Significance

dataForVolcanoPlot$Significance <- 'NS'

dataForVolcanoPlot$Signatures <- signatures$Common[match(dataForVolcanoPlot$Ensembl, signatures$`Ensembl ID`)]
dataForVolcanoPlot$Signature.Count <- signatures$Count[match(dataForVolcanoPlot$Ensembl, signatures$`Ensembl ID`)]


#saveRDS(dataForVolcanoPlot, file='manuscript/Figure_dataForVolcanoPlot_TCGA.RDS')

dataForVolcanoPlot <- readRDS(file='figures/Figure_dataForVolcanoPlot_TCGA.RDS')

cols <- c('UP'=google.red, 'NS'='gray','DOWN'=google.green)

dataForVolcanoPlot <- dataForVolcanoPlot[which(dataForVolcanoPlot$Signature.Count>=1),]

#View(dataForVolcanoPlot)



p <- ggplot(dataForVolcanoPlot, aes(x = logFC, y = -log10(adj.P.Val))) +
  #xlim(-2,2) +
  labs(x=expression(bold('Log'['2']*'(Fold Change)')), 
       y=(expression(bold('-Log'['10']*'(FDR)'))), 
       title=NULL) +
  geom_point(color='gray', alpha=1, size=2) +
  geom_vline(xintercept = hrThreshold,
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold), 
             color='darkgreen',linetype='dashed')+
  #xlim(-1,1) +
  #scale_x_continuous(breaks=c(0,1,2,4,6)) +
  #scale_y_continuous(expand = c(0.3, 0)) +
  #scale_color_manual(values = c('#4285F4',"gray", '#FBBC05')) +
  scale_color_manual(values = cols) +
  #facet_wrap(~Comparison, ncol = 2) +
  # geom_point(data = subset(dataForVolcanoPlot,
  #                          miRNA.ID %in% c('hsa-miR-7-1-3p','hsa-miR-7-5p')),
  #            aes(x=logFC, y=-log10(adj.P.Val)), size=2, color='blue', fill='gray',shape=21)+
  
  
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance.Label == 'UP' & Signature.Count>=1),
             aes(x=log2(HR), y=-log10(FDR)), fill=google.red, 
             size=2, color=google.red,shape=21)+
  
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance.Label == 'DOWN' & Signature.Count>=1),
             aes(x=log2(HR), y=-log10(FDR)), fill=google.green, 
             size=2, color=google.green,shape=21)+
  scale_fill_manual(values = cols) +
  
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance.Label == 'UP' & Signature.Count>=3),
                  aes(label = Symbol), 
                  nudge_x=0.25,
                  #force=1,
                  segment.alpha = 0.5, #segment.size = 0.5,
                  min.segment.length = 0,
                  size = 4, color='black', segment.color = 'black') +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance.Label == 'DOWN' & Signature.Count>=3),
                  aes(label = Symbol),
                  nudge_x=-0.2,
                  #force=1,
                  segment.alpha = 0.5, #segment.size = 0.5,
                  #min.segment.length = 5,
                  size = 4, color='black', segment.color = 'black') +
  theme_bw() +
  theme(axis.line = element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        strip.text = element_text(size=14, face='bold')) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))







p <- ggplot(dataForVolcanoPlot, aes(x = log2(HR), y = -log10(FDR))) +
  #xlim(-2,2) +
  labs(x=expression(bold('Log'['2']*'(Hazard Ratio)')), 
       y=(expression(bold('-Log'['10']*'(FDR)'))), 
       title=NULL) +
  geom_point(color='gray', alpha=1, size=2) +
  geom_vline(xintercept = hrThreshold,
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold), 
             color='darkgreen',linetype='dashed')+
  #xlim(-1,1) +
  #scale_x_continuous(breaks=c(0,1,2,4,6)) +
  #scale_y_continuous(expand = c(0.3, 0)) +
  #scale_color_manual(values = c('#4285F4',"gray", '#FBBC05')) +
  scale_color_manual(values = cols) +
  #facet_wrap(~Comparison, ncol = 2) +
  # geom_point(data = subset(dataForVolcanoPlot,
  #                          miRNA.ID %in% c('hsa-miR-7-1-3p','hsa-miR-7-5p')),
  #            aes(x=logFC, y=-log10(adj.P.Val)), size=2, color='blue', fill='gray',shape=21)+
  
  
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance.Label == 'UP' & Signature.Count>=1),
             aes(x=log2(HR), y=-log10(FDR)), fill=google.red, 
             size=2, color=google.red,shape=21)+
  
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance.Label == 'DOWN' & Signature.Count>=1),
             aes(x=log2(HR), y=-log10(FDR)), fill=google.blue, 
             size=2, color=google.blue,shape=21)+
  scale_fill_manual(values = cols) +
  
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance.Label == 'UP' & Signature.Count>=3),
                  aes(label = Symbol), 
                  nudge_x=0.1,
                  #force=1,
                  segment.alpha = 0.5, #segment.size = 0.5,
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  size = 4, color='black', segment.color = 'black') +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance.Label == 'DOWN' & Signature.Count>=3),
                  aes(label = Symbol),
                  nudge_x=-0.4,
                  #force=1,
                  segment.alpha = 0.5, #segment.size = 0.5,
                  max.overlaps = Inf,
                  #min.segment.length = 5,
                  size = 4, color='black', segment.color = 'black') +
  theme_bw() +
  theme(axis.line = element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        strip.text = element_text(size=14, face='bold')) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


p


topptx(p, filename = 'figures/Figure_Volcano_CoxPH_blue_red.pptx', width=6, height = 6)



##############################

max(dataForForestPlot$HR)

plotForest <- function(df, var.name="Drug", ratio.name="HR", ylab="", nlab="N (%)"){
  google.red <- '#ea4235'
  google.yellow <- '#fabd03'
  google.green <- '#34a853'
  google.blue <- '#4286f5'
  
  df[, ratio.name] <- round(df[, ratio.name], 2)
  df$Lower95 <- round(df$Lower95, 2)
  df$Upper95 <- round(df$Upper95, 2)
  df$P <- sapply(df$P, function(x) ifelse(x<0.001, formatC(x, format="e", digits=2), round(x, 3)))
  
  df$N <- paste0(df$Lower95, ' - ', df$Upper95)
  
  rangeb <- c(0,max(df$Upper95))
  rangeplot <- rangeb
  rangeplot[1] <- rangeb[2]*-1*1.3
  rangeplot[2] <- rangeplot[2] + 0.01 * diff(rangeb)
  cpositions=0.075*c(0.8, 2.5, 4.5, 6.5)
  width <- diff(rangeplot)
  # y-coordinates for labels:
  y_variable <- rangeplot[1] +  cpositions[1] * width
  y_nlevel <- rangeplot[1] +  cpositions[2] * width
  y_n <- rangeplot[1]  +  cpositions[3] * width
  y_cistring <- rangeplot[1]  +  cpositions[4] * width
  y_stars <- rangeb[2]
  plt_title <- ifelse(ratio.name=="OR", "Odds Ratio\n(95% CI)", "Hazard Ratio")
  
  p <- ggplot(df, aes_string(x=seq_along(df[, var.name]), y=ratio.name)) +
    geom_rect(aes(xmin=seq_along(df[, var.name])-0.5, xmax=seq_along(df[, var.name])+0.5,
                  ymin=rangeplot[1], ymax=rangeplot[2],
                  fill=ordered(seq_along(df[, var.name]) %% 2 + 1))) +
    scale_fill_manual(values = c("#00000033", "#FFFFFF33"), guide = "none") +
    scale_y_continuous(breaks = c(0,0.5,1,1.5,2), 
                       labels = c(0,0.5,1,1.5,2)) +
    scale_x_continuous(limits = c(0, length(df[,var.name])+2.2), expand = c(0,0)) +
    geom_errorbar(aes(ymin=Lower95, ymax=Upper95),width=0.4, size=0.8, color='black')+ 
    geom_point(color=google.red, size=3, shape=15) + 
    geom_text(data=df, aes(x=seq_along(df[, var.name]), y=y_variable, label=df[, var.name], group=NULL),
              size=5) +
    geom_text(data=df, aes(x=seq_along(df[, var.name]), y=y_nlevel, label=df[, ratio.name], group=NULL),
              size=5) +
    geom_text(data=df, aes(x=seq_along(df[, var.name]), y=y_n, label=df[, "N"], group=NULL),
              size=5) +

    # geom_text(data=df, aes(x=seq_along(df[, var.name]), y=y_nlevel, label=paste0('(', Lower95, '-', Upper95, ')'), group=NULL),
    #           size=5, vjust=1.5) +
    geom_text(data=df, aes(x=seq_along(df[, var.name]), y=y_cistring, label=P, group=NULL),
              size=5) +
    geom_hline(yintercept=1, linetype=2, color='black') +
    geom_hline(yintercept=0, linetype=1, color='grey') +
    #scale_x_reverse() + 
    coord_flip()+
    xlab('')+
    ylab(ylab) +
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.text = element_text(size=14),
          legend.position = 'right') +
    theme(axis.title.x=element_text(size=16, face = 'bold', hjust = 0.71),
          axis.title.y=element_blank(),
          axis.ticks=element_blank(),
          axis.text = element_text(color='black', size=14, face = 'bold'),
          axis.text.y = element_blank(),
          strip.text = element_text(size=14)) +
    theme(
      axis.line.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()) +
    annotate(geom='text', x=length(df[,var.name])+1.2, y=c(y_variable,y_nlevel, y_n, y_cistring),
             label=c('Gene', 'Hazard Ratio', '95% CI', 'P Value'), size=5, fontface='bold')
  return(p)
}



dataForForestPlot <- dataForVolcanoPlot
dataForForestPlot <- dataForForestPlot[which(dataForForestPlot$Signature.Count>=3),]
dim(dataForForestPlot)

dataForForestPlot <- dataForForestPlot[rev(rownames(dataForForestPlot)),]

#dataForForestPlot <- dataForForestPlot[order(dataForForestPlot$HR, decreasing = T),]


plotForest(dataForForestPlot, var.name="Symbol", ratio.name="HR", ylab="", nlab="Signature.Count")


########################################################


library(clusterProfiler)
library(DOSE)
library(ReactomePA)

enrichmentFun <- function(targets, ont, gene.ids) {
  
  if (ont=='DO') {
    
    kk <- enrichDO(gene          = targets,
                   ont           = 'DO',
                   pAdjustMethod = 'BH',
                   pvalueCutoff  = 1,
                   minGSSize     = 10,
                   maxGSSize     = 500,
                   readable      = FALSE)
    
    
  } else if (ont=='KEGG') {
    
    kk <- enrichKEGG(gene          = targets,
                     organism      = 'hsa',
                     pAdjustMethod = 'BH',
                     pvalueCutoff  = 1,
                     minGSSize     = 10,
                     maxGSSize     = 500)
    
    
  } else if (ont == 'HALLMARK') {
    
    read.gmt = function(file){
      if(!grepl("\\.gmt$",file)[1]){stop("Pathway information must be a .gmt file")}
      geneSetDB = readLines(file)                                ##read in the gmt file as a vector of lines
      geneSetDB = strsplit(geneSetDB,"\t")                       ##convert from vector of strings to a list
      names(geneSetDB) = sapply(geneSetDB,"[",1)                 ##move the names column as the names of the list
      geneSetDB = lapply(geneSetDB, "[",-1:-2)                   ##remove name and description columns
      geneSetDB = lapply(geneSetDB, function(x){x[which(x!="")]})##remove empty strings
      return(geneSetDB)
    }
    
    gmtfile <- 'manuscript/h.all.v7.3.entrez.gmt'
    h <- read.gmt(gmtfile)
    
    h <- data.frame(gs_name=rep(names(h), as.numeric(lapply(h, length))),
                    entrez_gene=as.character(unlist(h)))
    
    kk <- enricher(gene          = targets,
                   TERM2GENE     = h,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1,
                   minGSSize     = 10,
                   maxGSSize     = 500)
    
  } else if (ont=='REACTOME') {
    
    kk <- enrichPathway(gene          = targets,
                        organism      = 'human',
                        pAdjustMethod = 'BH',
                        pvalueCutoff  = 1,
                        minGSSize     = 10,
                        maxGSSize     = 500,
                        readable      = FALSE)
    
  } else if (ont=='GOBP') {
    
    kk <- enrichGO(gene          = targets,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1,
                   minGSSize     = 10,
                   maxGSSize     = 500,
                   readable      = FALSE)

  }
  
  kk <- kk@result
  
  List.Total <- length(targets)
  Count <- kk$Count
  Pop.Hits <- as.numeric(sapply(kk$BgRatio, function(x) strsplit(x, '/')[[1]][1]))
  Pop.Total <- as.numeric(sapply(kk$BgRatio, function(x) strsplit(x, '/')[[1]][2]))
  
  Fold.Enrichment <- format(Count/List.Total*Pop.Total/Pop.Hits, digits=2, nsmall=2)
  
  P.Value <- format(kk$pvalue, digits=3, nsmall=3)
  BH.Adj.P <- formatC(kk$p.adjust, format = 'e', digits=2)
  
  Entrez <- kk$geneID
  genes <- strsplit(kk$geneID, '/')
  
  
  #gene.ids <- signatures
  
  idx <- lapply(genes, function(v) match(v, gene.ids$`Entrez ID`))
  Ensembl <- unlist(lapply(idx, function(x) paste(gene.ids$`Ensembl ID`[x], collapse = '/')))
  Symbol <- unlist(lapply(idx, function(x) paste(gene.ids$`HGNC Symbol`[x], collapse = '/')))
  
  enrichTable <- data.frame(ID = kk$ID,
                   Description = kk$Description,
                   Count, List.Total, Pop.Hits, Pop.Total,
                   Fold.Enrichment, P.Value, BH.Adj.P,
                   Entrez, Ensembl, Symbol,
                   stringsAsFactors = F)
  
  
  return (enrichTable)
  
  
}


idx <- match(gene.stats$`Ensembl ID`[gene.stats$Count>=3], signatures$`Ensembl ID`)
signatures$`Entrez ID`[idx]
targets <- signatures$`Entrez ID`[idx]

enrichTable <- enrichmentFun(targets = targets, ont = 'DO', gene.ids = signatures)

saveRDS(enrichTable, file='manuscript/GOBP_40_Genes.RDS')


###
idx <- match(gene.stats$`Ensembl ID`[gene.stats$Count>=1], signatures$`Ensembl ID`)
signatures$`Entrez ID`[idx]
targets <- signatures$`Entrez ID`[idx]


# 40
# 142
# 1032

for (ont in c('DO','KEGG','HALLMARK','REACTOME','GOBP')) {
  
  message (ont)
  
  enrichTable <- enrichmentFun(targets = targets, ont = ont, gene.ids = signatures)
  
  saveRDS(enrichTable, file=paste0('report/Enrichment/', ont, '_142_Genes.RDS'))
  
}

unique(signatures$Signature)




for (ont in c('DO','KEGG','HALLMARK','REACTOME','GOBP')) {
  
  message (ont)
  
  enrichList <- list()
  
  
  for (.sig in unique(signatures$Signature)) {
  
    print (.sig)
      
    .idx <- which(signatures$Signature==.sig)
    
    targets <- signatures$`Entrez ID`[.idx]
    
    enrichTable <- enrichmentFun(targets = targets, ont = ont, gene.ids = signatures)
  
    enrichList[[.sig]] <- enrichTable
    
  }
  
  saveRDS(enrichList, file=paste0('report/Enrichment/', ont, '_Signature_Genes.RDS'))
  
}



EnrichmentBubblePlotFun <- function(dataForBubblePlot, topn = 30) {
  
  if (nrow(dataForBubblePlot)>topn) {
    dataForBubblePlot <- dataForBubblePlot[1:topn,]
  }
  
  p <- ggplot(dataForBubblePlot, mapping=aes(x=Description, y=Fold.Enrichment, #y=-log10(Benjamini), #y=Fold.Enrichment
                                             color=BH.Adj.P,size=Count)) +
    geom_point()+ coord_flip() +
    scale_x_discrete(limits=rev(unique(dataForBubblePlot$Description))) +
    #scale_x_discrete(limits=Order)+
    scale_colour_gradientn(limits=c(0,0.05),
                           colors= c("red","yellow","green")) + #
    #facet_wrap(~Comparison) +
    #facet_grid(Regulation~Comparison) + # scales=free
    xlab('')+ylab('Fold Enrichment') + #ggtitle("") + 
    guides(shape = guide_legend(order=1),
           colour = guide_colourbar(order=2, title = 'FDR')) + #'P Value\n(Benjamini)'))
    theme_bw()+theme(axis.line = element_line(colour = "black"),
                     panel.grid.minor = element_blank(),
                     panel.border = element_rect(colour='black'),
                     panel.background = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5, size=20)) +
    theme(axis.text=element_text(size=14, color='black', face = 'bold'),
          axis.text.x =element_text(size=14, color='black', face = 'bold', angle=0, hjust=0.5),
          axis.title=element_text(size=16, face = 'bold')) +
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    theme(#strip.text = element_text(size = 14),
      legend.key.size = unit(0.8,'cm'))
  
  
  return (p)
  
}

dataForBubblePlot <- readRDS('report/Enrichment/KEGG_142_Genes.RDS')


enrichList <- readRDS('report/Enrichment/KEGG_Signature_Genes.RDS')

names(signature.list)

dataForBubblePlot <- enrichList[['Wu']]

dataForBubblePlot <- readRDS('figures/DO_40_Genes.RDS')

dataForBubblePlot

dataForBubblePlot$Fold.Enrichment <- as.numeric(dataForBubblePlot$Fold.Enrichment)
dataForBubblePlot$Count <- as.numeric(dataForBubblePlot$Count)
dataForBubblePlot$BH.Adj.P <- as.numeric(dataForBubblePlot$BH.Adj.P)

#dataForBubblePlot <- dataForBubblePlot[dataForBubblePlot$BH.Adj.P<=0.05,]

capitalize <- function(x) {
  substr(x,1,1) <- toupper(substr(x,1,1))
  return(x)
}

dataForBubblePlot$Description <- capitalize(dataForBubblePlot$Description)
dataForBubblePlot <- dataForBubblePlot[which(dataForBubblePlot$BH.Adj.P<0.05),]

#View(dataForBubblePlot)

p <- EnrichmentBubblePlotFun(dataForBubblePlot)

topptx(p, filename = 'figures/Figure_DO_40_Genes.pptx', width=8.5, height = 6)


################################################

.g <- gene.stats$`Ensembl ID`[gene.stats$Count>=3]
.g


dim(exprData)

expr <- exprData[.g,]
rownames(expr) <- gene.stats$`HGNC Symbol`[gene.stats$Count>=3]

M <- cor(t(expr))
M

library(corrplot)

col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))

corrplot(M, type = "upper", 
         order = c("original"),
         col = rev(col2(200)),
         tl.col = "black",
         tl.cex = 1, 
         #mar=c(0.1,0.1,0.1,0.1)
         )

# corrplot.mixed(M, #type = "upper", 
#          order = c("original"),
#          number.cex = .5,
#          
#          #col = rev(col2(200))#,
#          #tl.col = "black",
#          tl.cex = 0.5
#          #mar=c(0.1,0.1,0.1,0.1)
# )





