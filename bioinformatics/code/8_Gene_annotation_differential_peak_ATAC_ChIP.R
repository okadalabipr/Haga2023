# 8_Gene_annotation_differential_peak_ATAC_ChIP
# Gene annotation using ChIPseeker
# Generate Venn figure of RNA-seq, ATAc-seq, and ChIP-seq (H3K27Ac)
# Generate Venn figure of RNA-seq(Upregulated) and ChIP-seq (H3K27Ac gained)

# Libraies I need
library(tidyverse)
library(clusterProfiler)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(imguR)
library(UpSetR)
library(VennDiagram)
library(ggplotify)

#Set were to save data
setwd("D:/haga/HFF_RS_TGF_ChIP/R/data/ChIPseeker/Paper/") 

# Check were you are now
getwd() 

# clear the decks
rm(list = ls())

# load hg38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# H3K27Ac
# read Differential peak files(use nf-core output)
H3K27Ac_p_0.05 <- readPeakFile("D:/haga/HFF_RS_TGF_ChIP/results/bwa/mergedLibrary/macs/broadPeak/consensus/Diff_H3K27Ac_PDL/Differential_H3K27Ac_output.bed")
# ATAC
# read Differential peak files(use nf-core output)
ATAC_p_0.05 <- readPeakFile("D:/haga/HFF_RS_TGF_ATAC/results/bwa/mergedLibrary/macs/broadPeak/consensus/Diff_ATAC_PDL/Differential_ATAC_output.bed")

#list
# Peak Annotation
peakAnno_H3K27Ac_p_0.05 <- annotatePeak(H3K27Ac_p_0.05,
                                        tssRegion=c(-3000, 3000),
                                        TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_ATAC_p_0.05 <- annotatePeak(ATAC_p_0.05,
                                     tssRegion=c(-3000, 3000),
                                     TxDb=txdb, annoDb="org.Hs.eg.db")

#Obtain ENSEMBL gene ID
# ChIp-seq(H3K27Ac)
#1185 genes
ChIPseq <- peakAnno_H3K27Ac_p_0.05@anno$ENSEMBL %>% as.data.frame() %>% distinct() %>% dplyr::rename("ENSEMBL" = ".") %>% na.omit()
sum(table(ChIPseq))
ChIPseq_ID <- ChIPseq$ENSEMBL
length(ChIPseq_ID)

#ATAC-seq
#1254 genes
ATACseq <- peakAnno_ATAC_p_0.05@anno$ENSEMBL %>% as.data.frame() %>% distinct() %>% dplyr::rename("ENSEMBL" = ".") %>% na.omit()
sum(table(ATACseq))
ATACseq_ID <- ATACseq$ENSEMBL
length(ATACseq_ID)

#Import RNA-seq DEGs
PDL24_vs_PDL36_RNA <- read.csv("D:/haga/HFF_RS_TGF/R/data_R/150bp/DESeq2/PDL24_vs_PDL36_DEGs_DoRothEA.csv",header = T,stringsAsFactors = FALSE)
PDL24_vs_PDL47_RNA <- read.csv("D:/haga/HFF_RS_TGF/R/data_R/150bp/DESeq2/PDL24_vs_PDL47_DEGs_DoRothEA.csv",header = T,stringsAsFactors = FALSE)
PDL36_vs_PDL47_RNA <- read.csv("D:/haga/HFF_RS_TGF/R/data_R/150bp/DESeq2/PDL36_vs_PDL47_DEGs_DoRothEA.csv",header = T,stringsAsFactors = FALSE)

#For RNAseq, integrate all DEGs in RNAseq_ID
RNAseq <- full_join(PDL24_vs_PDL36_RNA,
                    PDL24_vs_PDL47_RNA,
                    by="ensembl_gene_id") 
RNAseq <- full_join(RNAseq,
                    PDL36_vs_PDL47_RNA,by="ensembl_gene_id") 
RNAseq_ID <- RNAseq$ensembl_gene_id
length(RNAseq_ID)

#Venn
color_palette <- c("#4F81BD","#C00000","#BC80BD")
list <- list(RNAseq =RNAseq_ID,
             H3K27Ac =ChIPseq_ID,
             ATACseq =ATACseq_ID )
venn.diagram(list,
             filename="Fig1_G_Venn_intersect_RNA_H3K27Ac_ATAC_0.05.png",
             category.names = c("RNA-seq","H3K27Ac\nChIP-seq","ATAC-seq"),
             imagetype="png" ,
             height = 480 , 
             width = 480 , 
             resolution = 300,
             lwd = 2,
             lty = 'blank',
             fill = color_palette,
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             cat.cex = 0.6,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             alpha=0.4)

#Get intersected gene list
genes_merge <- intersect(RNAseq_ID,ChIPseq_ID)
genes_merge <- intersect(genes_merge,ATACseq_ID)
genes_merge_dataframe <-genes_merge %>% as.data.frame() %>% distinct() %>% dplyr::rename("ENSEMBL" = ".") %>% na.omit()

#Save genes for IPA analysis
write.csv(genes_merge_dataframe, file = "RNAseq_H3K27Ac_ATACseq_merge_genename.csv")

#------------------------------------------------------------------------------------------------------------
# Code for Figure S2A
# H3K27Ac
# read peak files for gained H3K27Ac peak
Gain_H3K27Ac_p_0.05 <- readPeakFile("D:/haga/HFF_RS_TGF_ChIP/results/bwa/mergedLibrary/macs/broadPeak/consensus/H3K27Ac/deseq2/H3K27Ac_output_gain.bed")

#list
# Peak Annotation
peakAnno_Gain_H3K27Ac_p_0.05 <- annotatePeak(Gain_H3K27Ac_p_0.05,
                                             tssRegion=c(-3000, 3000),
                                             TxDb=txdb, annoDb="org.Hs.eg.db")

#Obtain ENSEMBL gene ID
peakAnno_Gain_H3K27Ac_p_0.05_GENE <- peakAnno_Gain_H3K27Ac_p_0.05@anno$ENSEMBL %>% as.data.frame() %>% distinct() %>% dplyr::rename("ENSEMBL" = ".") %>% na.omit()

#Import RNA-seq DEGs
upregulated_RNA <- read.csv("D:/haga/HFF_RS_TGF/R/data_R/150bp/DESeq2/20221008_HFF_RS_DEG_up.csv",header = T,stringsAsFactors = FALSE)

#For RNAseq
#upregulated genes 1155
RNAseq_up_ID <- upregulated_RNA$ENSEMBL
length(RNAseq_up_ID)

#H3K27Ac
#Gain gene 1093
H3K27Ac_gain_ID <- peakAnno_Gain_H3K27Ac_p_0.05_GENE$ENSEMBL
length(H3K27Ac_gain_ID)

#Venn
color_palette <- c("#4F81BD","#C00000")
list <- list(RNAseq =RNAseq_up_ID,
             H3K27Ac =H3K27Ac_gain_ID)

venn.diagram(list,
             filename="FigS2_A_Venn_intersect_Upregulate_H3K27Ac_gain.png",
             category.names = c("RNAseq" ,"H3K27Ac"),
             imagetype="png" ,
             height = 480 , 
             width = 480 , 
             resolution = 300,
             lwd = 2,
             lty = 'blank',
             fill = color_palette,
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             cat.cex = 0.6,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27),
             cat.dist = c(0.055, 0.055),
             cat.fontfamily = "sans",
             alpha=0.4)

#Get gene list
genes_merge <- intersect(RNAseq_up_ID,H3K27Ac_gain_ID)
genes_merge <- intersect(genes_merge,ATAC_gain_ID)
genes_merge_dataframe <-genes_merge %>% as.data.frame() %>% distinct() %>% dplyr::rename("ENSEMBL" = ".") %>% na.omit()

#Save
#Use for IPA analysis
write.csv(genes_merge_dataframe, file = "Upregulated_RNAseq_H3K27Ac_gain_merge_genename.csv")
