# 11_RNA_ATAC_Con_vs_TGF_THBS1_FMOD
# Analysis of THBS1 and FMOD using RNA-seq and ATAC-seq

# Libraies I need
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(ggrepel)
library(export)

# > packageVersion("ChIPseeker")
# [1] 1.32.1
packageVersion("ChIPseeker")
# > packageVersion("ClusterProfiler")
# [1] 4.9.1
packageVersion("ClusterProfiler")

#Set were to save data
setwd("D:/haga/github/Haga2023/bioinformatics/ref_file/") 

# Check were you are now
getwd() 

# clear the decks
rm(list = ls())

# Make object
counts_deseq <- read.table("20211025_HFF_counts.txt", header=T, stringsAsFactors = F)

# remove unwanted lines
counts_deseq <- counts_deseq[,!colnames(counts_deseq) %in% c("Chr", "Start", "End", "Strand")]
counts_deseq
counts_deseq <- counts_deseq[,c(1:8)]
colnames(counts_deseq)
colnames(counts_deseq) <- c("Geneid","Length",
                            "HFF_PDL24_CON_1", "HFF_PDL24_CON_2", "HFF_PDL24_CON_3",
                            "HFF_PDL24_TGF_1", "HFF_PDL24_TGF_2", "HFF_PDL24_TGF_3")
#ion
counts_deseq[,1] <- str_replace(counts_deseq[,1],
                                pattern = ".[0-9]+$",
                                replacement = "")

counts_deseq <- counts_deseq[,c(1,3:ncol(counts_deseq))]
colnames(counts_deseq)
head(counts_deseq)

# import data
gene_set <- read.csv("TPM_normalized_data_cutoff_mean_5_proteincoding.csv", header = T, stringsAsFactors = F)
gene_set <- gene_set$X %>% as.data.frame()
gene_set <- dplyr::rename(gene_set, Geneid = .)
gene_set

# gene_set(TPM>5)
#10303 genes
sum(table(gene_set))
counts_deseq <- dplyr::inner_join(gene_set, counts_deseq, by = "Geneid")
colnames(counts_deseq)
dim(counts_deseq)

#group
Con_TGF <-counts_deseq %>% column_to_rownames(var = "Geneid")
colnames(Con_TGF)
#------------------------------------------------------------------------------------------------------------

#Make Tag for each groups
group_Con_TGF <- data.frame(con = factor(c("Con","Con","Con",
                                           "TGF","TGF","TGF")))

# Change row name
rownames(group_Con_TGF) <- c(colnames(Con_TGF))

# Make_data_set_for_deseq
dds_1 <- DESeqDataSetFromMatrix(countData = Con_TGF, colData = group_Con_TGF, design = ~ con)

# Set ref
dds_1$con <- relevel(dds_1$con, ref = "Con")

# check
relevel(dds_1$con, ref = "Con")

# Run_deseq
dds_1 <- DESeq(dds_1)

# result
res_A <- results(dds_1)

# Make dataset for p value
res_A<- res_A[,c(2,6)]%>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")

#rename
colnames(res_A) <- c("ENSEMBL","log2FoldChange","padj")

# identify data of p<0.05
res_A$sig <- if_else(res_A$padj <0.05, "TRUE", "FALSE")
table(res_A$sig)

# padj<0.05, 
#arrage by log2FC
graph_1 <- res_A %>% dplyr::arrange(by = log2FoldChange)
graph_1
#Save DEGs
write.csv(graph_1,file = "PDL24_Con_vs_TGF_DEGs_p_0.05.csv")

# load hg38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#Read RNA-seq
#Read Deseq2 outputs
#PDL24
#Control vs TGF-beta1
PDL24_RNAseq <- read.csv("PDL24_Con_vs_TGF_DEGs_p_0.05.csv",header = T,row.names = 1)
PDL24_RNAseq <- PDL24_RNAseq %>% filter(padj<0.05)
PDL24_RNAseq <- PDL24_RNAseq[,c(1,2)]

#6118 genes from RNA-seq
nrow((PDL24_RNAseq))

#ATAC
#Read Deseq2 outputs
#PDL24
#Control vs TGF-beta1
#FDR < 0.05
PDL24 <- read.table("HFF_PDL_24_CONvsHFF_PDL_24_TGF.mLb.clN.deseq2.FDR0.05.results.txt",header = T)
PDL24 <- PDL24[,c(2:6,8,12)]
PDL24 <- PDL24 %>% filter(padj<0.05)
head(PDL24)

#Save
write.table(PDL24,"ATAC_PDL24_Control_vs_TGF_padj_0.05.bed",quote=F,sep="\t",row.names=F,col.names=F)

# ATAC
# read peak files
ATAC_PDL24_p_0.05 <- readPeakFile("ATAC_PDL24_Control_vs_TGF_padj_0.05.bed")
ATAC_PDL24_p_0.05
#list
# Peak Annotation
peakAnno_ATAC_PDL24_p_0.05 <- annotatePeak(ATAC_PDL24_p_0.05,
                                           tssRegion=c(-3000, 3000),
                                           TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_ATAC_PDL24_p_0.05_anno <- peakAnno_ATAC_PDL24_p_0.05@anno  %>% as.data.frame()
colnames(peakAnno_ATAC_PDL24_p_0.05_anno)
PDL24_ATAC <- peakAnno_ATAC_PDL24_p_0.05_anno[,c(19,8)] %>% na.omit()
PDL24_ATAC

#Join RNA-seq and ATAC-seq
Join_RNAseq_ATACseq <- inner_join(PDL24_RNAseq,PDL24_ATAC,by="ENSEMBL") 
colnames(Join_RNAseq_ATACseq) <-c("ENSEMBL","RNAseq_log2FC","ATACseq_log2FC")
Join_RNAseq_ATACseq
PDL24_RNAseq

#Chage ENSEMBL to SYMBOL
ENSEMBL <- Join_RNAseq_ATACseq$ENSEMBL
SYMBOL <- bitr(ENSEMBL,  fromType = "ENSEMBL",
               toType = c("SYMBOL"),
               OrgDb = org.Hs.eg.db)

res <- inner_join(SYMBOL,Join_RNAseq_ATACseq,by="ENSEMBL") 
res

#label
res$THBS1_FMOD <- if_else(res$SYMBOL=="THBS1"|res$SYMBOL=="FMOD", "TRUE", "FALSE")
table(res$THBS1_FMOD)

#label to show
res_slice <- res %>% dplyr::filter(THBS1_FMOD=="TRUE")
selected_gene <- res_slice$SYMBOL
labels_to_show <- res$SYMBOL
labels_to_show[!labels_to_show %in% selected_gene] <- ""
labels_to_show

#Plot
color_palette <- c("gray50","#C00000")
g <- ggplot(res, aes(x = RNAseq_log2FC, y = ATACseq_log2FC, label =labels_to_show, color = THBS1_FMOD)) 
g <- g + geom_point(size =0.8, alpha =0.5) 
g <- g + theme_bw() 
g <- g + scale_color_manual(values = color_palette)
g <- g + xlab (expression("RNA-seq (log2FC)")) + ylab(expression("ATAC-seq (log2FC)")) 
g <- g + geom_text_repel(size=2,max.overlaps = getOption("ggrepel.max.overlaps", default = 1000)) 
g <- g + guides(colour=FALSE) 
g <- g + theme(axis.title = element_text(size = 7), 
               axis.text = element_text(size = 7), 
               axis.text.x = element_text(size = 7), 
               axis.text.y = element_text(size = 7), 
               legend.text = element_text(size = 7), 
               legend.title = element_text(size = 7)) + theme(plot.title = element_text(size = 7))
g

graph2ppt(x=g, file="Fig1_J_RNAseq_ATACseq_LOg2FC_Plot", width = 3, height = 2.2)
