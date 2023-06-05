# 4_Transcription_factor_enrichment
# DoRothEA
# Identfication of TFs by RNA-seq

#Tutorial used for analysis
#https://github.com/saezlab/transcriptutorial/blob/master/scripts/04_TranscriptionFactor_activity_with_Dorothea.md

# Libraies I need
library("progeny")
library("dorothea")
library("tidyverse")
library("pheatmap")
library("readr")
library("ggrepel")
library("tidyverse")
library("data.table")
library("edgeR")
library("biomaRt")
library("ggrepel")
library("ggsci")
library("DESeq2")
library("gdata") 

# Set were to save data
setwd("d:/haga/HFF_RS_TGF/R/data_R/150bp/DoRothEA/") 

## We also load the support functions
source("support_functions.R")

# Check were you are now
getwd() 

# clear the decks
rm(list = ls())

#Calculate TPM value
#Read output from feautureCounts
#Make object d
d <- read.table("d:/haga/HFF_RS_TGF/featureCounts/20211025_HFF_counts.txt", header=T, stringsAsFactors = F)

# remove unwanted lines
d <- d[,!colnames(d) %in% c("Chr", "Start", "End", "Strand")]

# Replicative stress
# PDL only
d <- d[c(1:5,9:11,15:17)]
colnames(d)
# change names of samples(initial names are long)
colnames(d) <- c("Geneid","Length",
                            "HFF_PDL24_CON_1", "HFF_PDL24_CON_2", "HFF_PDL24_CON_3",
                            "HFF_PDL36_CON_1", "HFF_PDL36_CON_2", "HFF_PDL36_CON_3",
                            "HFF_PDL47_CON_1", "HFF_PDL47_CON_2", "HFF_PDL47_CON_3")

# Divide the sample read count by the length of the gene multiplied by 1000
d_tpm <- (d[,3:ncol(d)]/d$Length) *1000
d_tpm
# combind Geneid from d to d_tpm
d_tpm <- cbind(ensembl_gene_id = d$Geneid, d_tpm)

head(d_tpm)

# Get the normalization factor, if group doesn't exist, will be in same group
normfactor<- DGEList(counts = d_tpm[,2:ncol(d_tpm)], group = colnames(d_tpm[,2:ncol(d_tpm)]))
normfactor <- calcNormFactors(normfactor, method="RLE")
normfactor$counts
normfactor$samples

# - $counts 
normfactor_samples <- normfactor$samples
normfactor_samples <- as.data.frame(normfactor_samples)
normfactor_samples

#write graph
g <- ggplot(normfactor_samples, aes(x =group, y = lib.size)) 
g <- g + geom_point(size =4, alpha =0.5) + theme_bw() + xlab("sample")+ ylab("lib.size") + theme(axis.title = element_text(size = 16))
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))
g <- g + theme(axis.title = element_text(size = 18), 
               axis.text.x = element_text(size = 14), 
               axis.text.y = element_text(size = 14), 
               plot.title = element_text(size = 18))
g <- g + scale_y_continuous(expand = c(0, 0), limits = c(0,6000000))
g

# multiply normalization factor with the library size
normfactor_samples$normlib <- normfactor_samples$lib.size*normfactor_samples$norm.factors

# Loop to get the final TPM counts table
for(i in 1:(dim(d_tpm)[2]-1)){
  
  d_tpm[,i+1] <- (d_tpm[,i+1]/normfactor_samples$normlib[i])*1000000
  
}

d_tpm

# Check effect of TPM normalization
normfactor_after<- DGEList(counts = d_tpm[,2:ncol(d_tpm)], group = colnames(d_tpm[,2:ncol(d_tpm)]))
normfactor_after <- calcNormFactors(normfactor_after, method="RLE")

# - $counts 
normfactor_samples_after <- normfactor_after$samples

#write graph
g <- ggplot(normfactor_samples_after, aes(x =lib.size, y = norm.factors, label=group)) + geom_point(size =4, alpha =0.5) + theme_bw() + xlab("lib.size")+ ylab("norm.factors") + theme(axis.title = element_text(size = 16))
g <- g + geom_text_repel()
g

# Library read
# Use mmusculus dataset
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Show list of attributes 
listAttributes(mart)

#cut off version information
d_tpm[,1] <- str_replace(d_tpm[,1],
                         pattern = ".[0-9]+$",
                         replacement = "")
d_tpm

# Chage names of genes (add names of genes)
gene_list <-getBM(attributes = c("ensembl_gene_id","external_gene_name","gene_biotype"),
                  filters = 'ensembl_gene_id', 
                  values = d_tpm,
                  mart = mart)

# Combing gene_list & d_tpm
Final_gene_list <- dplyr::inner_join(gene_list, d_tpm ,by="ensembl_gene_id")
Final_gene_list
sum(table(Final_gene_list$ensembl_gene_id))

#Obtain genes TPM>5
norm_count_mean <- mutate (Final_gene_list, mean = rowMeans(Final_gene_list[,4:ncol(Final_gene_list)]))
norm_count_mean$sig <- if_else(norm_count_mean$mean>5, "TRUE", "FALSE")
norm_count_mean
# FALSE  TRUE 
# 48872 11533
table(norm_count_mean$sig)
sum(table(norm_count_mean$sig))
table(norm_count_mean$gene_biotype)

# Filter True
norm_count_mean_TRUE <- norm_count_mean %>% dplyr::filter(sig == TRUE)

#Protein conding 10148
table(norm_count_mean_TRUE$gene_biotype)
norm_count <- dplyr::select(norm_count_mean_TRUE, -sig & -mean)
sum(table(norm_count$external_gene_name))

# column to row names
norm_count <- column_to_rownames(norm_count, var = "ensembl_gene_id")

# Only protein coding genes
norm_count <- dplyr::filter(norm_count, gene_biotype == "protein_coding")

#Check how many genes
# protein_coding 
# 10303
table(norm_count$gene_biotype)

# save csv
write.csv(norm_count, file= "TPM_normalized_data_cutoff_mean_5_proteincoding.csv")
#------------------------------------------------------------------------------------------------------------
#Get DEGs between each PDL
#Again read count data
counts_deseq <- read.table("d:/haga/HFF_RS_TGF/featureCounts/20211025_HFF_counts.txt", header=T, stringsAsFactors = F)

# remove unwanted lines
counts_deseq <- counts_deseq[,!colnames(counts_deseq) %in% c("Chr", "Start", "End", "Strand")]

# PDL only
counts_deseq <- counts_deseq[c(1:5,9:11,15:17)]
colnames(counts_deseq)
# change names of samples(initial names are long)
colnames(counts_deseq) <- c("Geneid","Length",
                            "HFF_PDL24_CON_1", "HFF_PDL24_CON_2", "HFF_PDL24_CON_3",
                            "HFF_PDL36_CON_1", "HFF_PDL36_CON_2", "HFF_PDL36_CON_3",
                            "HFF_PDL47_CON_1", "HFF_PDL47_CON_2", "HFF_PDL47_CON_3")


#cut off version information
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
sum(table(gene_set))
counts_deseq <- dplyr::inner_join(gene_set, counts_deseq, by = "Geneid")
colnames(counts_deseq)
dim(counts_deseq)

#group
PDL24_PDL36 <-counts_deseq[,c(1,2:7)]%>% column_to_rownames(var = "Geneid")
PDL24_PDL47 <-counts_deseq[,c(1,2:4,8:10)]%>% column_to_rownames(var = "Geneid")
PDL36_PDL47 <-counts_deseq[,c(1,5:10)]%>% column_to_rownames(var = "Geneid")
colnames(PDL24_PDL36)
colnames(PDL24_PDL47)
colnames(PDL36_PDL47)

#Make Tag for each groups
group_PDL24_PDL36 <- data.frame(con = factor(c("PDL24","PDL24","PDL24","PDL36","PDL36","PDL36")))
group_PDL24_PDL47 <- data.frame(con = factor(c("PDL24","PDL24","PDL24","PDL47","PDL47","PDL47")))
group_PDL36_PDL47 <- data.frame(con = factor(c("PDL36","PDL36","PDL36","PDL47","PDL47","PDL47")))

# Change row name
rownames(group_PDL24_PDL36) <- c(colnames(PDL24_PDL36))
rownames(group_PDL24_PDL47) <- c(colnames(PDL24_PDL47))
rownames(group_PDL36_PDL47) <- c(colnames(PDL36_PDL47))

# Make_data_set_for_deseq
dds_1 <- DESeqDataSetFromMatrix(countData = PDL24_PDL36, colData = group_PDL24_PDL36, design = ~ con)
dds_2 <- DESeqDataSetFromMatrix(countData = PDL24_PDL47, colData = group_PDL24_PDL47, design = ~ con)
dds_3 <- DESeqDataSetFromMatrix(countData = PDL36_PDL47, colData = group_PDL36_PDL47, design = ~ con)

# Set ref
dds_1$con <- relevel(dds_1$con, ref = "PDL24")
dds_2$con <- relevel(dds_2$con, ref = "PDL24")
dds_3$con <- relevel(dds_3$con, ref = "PDL36")

# check
relevel(dds_1$con, ref = "PDL24")
relevel(dds_2$con, ref = "PDL24")
relevel(dds_3$con, ref = "PDL36")

# Run_deseq
dds_1 <- DESeq(dds_1)
dds_2 <- DESeq(dds_2)
dds_3 <- DESeq(dds_3)

# result
res_A <- results(dds_1)
res_B <- results(dds_2)
res_C <- results(dds_3)

# Make dataset for p value
res_A<- res_A[,c(2, 5)]%>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")
res_B<- res_B[,c(2, 5)]%>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")
res_C<- res_C[,c(2, 5)]%>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")

#merge files into one
res <- inner_join(res_A,res_B, by = "ENSEMBL")
res <- inner_join(res,res_C, by = "ENSEMBL")
colnames(res)
colnames(res) <- c("ENSEMBL",
                   "log2FoldChange_A","pvalue_A",
                   "log2FoldChange_B","pvalue_B",
                   "log2FoldChange_C","pvalue_C")
head(res)

# adjust p value among samples
res_A_p <- res$pvalue_A
res_B_p <- res$pvalue_B
res_C_p <- res$pvalue_C
length(res_A_p)

# Caculate FDR among samples
q <- p.adjust(c(res_A_p,res_B_p,res_C_p), method = "BH") %>% as.data.frame()
dim(q)

# extract each FDR
q_A <- q %>% dplyr::slice(1:10303) %>% as.data.frame()
q_B <- q %>% dplyr::slice(10304:20606) %>% as.data.frame()
q_C <- q %>% dplyr::slice(20607:nrow(q)) %>% as.data.frame()

#Check each gene number are same
dim(q_A)
dim(q_B)
dim(q_C)

#Integrate into one
res_A <- cbind(res$ENSEMBL,res$log2FoldChange_A,q_A) %>% dplyr::rename("padj" = ".")
res_B <- cbind(res$ENSEMBL,res$log2FoldChange_B,q_B) %>% dplyr::rename("padj" = ".")
res_C <- cbind(res$ENSEMBL,res$log2FoldChange_C,q_C) %>% dplyr::rename("padj" = ".")

#rename
colnames(res_A) <- c("ENSEMBL","log2FoldChange","padj")
colnames(res_B) <- c("ENSEMBL","log2FoldChange","padj")
colnames(res_C) <- c("ENSEMBL","log2FoldChange","padj")

# identify d#merge and outpu all data
temp_res <- dplyr::full_join(res_A,res_B) %>% full_join(res_C)
temp_res

#Data of p<0.01
res_A$sig <- if_else((res_A$log2FoldChange > log2(1.2) & res_A$padj <0.05) | (res_A$log2FoldChange < -log2(1.2) & res_A$padj <0.05), "TRUE", "FALSE")
res_B$sig <- if_else((res_B$log2FoldChange > log2(1.2) & res_B$padj <0.05) | (res_B$log2FoldChange < -log2(1.2) & res_B$padj <0.05), "TRUE", "FALSE")
res_C$sig <- if_else((res_C$log2FoldChange > log2(1.2) & res_C$padj <0.05) | (res_C$log2FoldChange < -log2(1.2) & res_C$padj <0.05), "TRUE", "FALSE")
table(res_A$sig)
table(res_B$sig)
table(res_C$sig)

res_B

#arrage by log2FC
graph_1 <- res_A %>% filter(sig == "TRUE")%>% dplyr::arrange(by = log2FoldChange)
graph_2 <- res_B %>% filter(sig == "TRUE")%>% dplyr::arrange(by = log2FoldChange)
graph_3 <- res_C %>% filter(sig == "TRUE")%>% dplyr::arrange(by = log2FoldChange)
graph_1

#Save DEGs
write.csv(graph_1,file = "PDL24_vs_PDL36_DEGs_DoRothEA.csv")
write.csv(graph_2,file = "PDL24_vs_PDL47_DEGs_DoRothEA.csv")
write.csv(graph_3,file = "PDL36_vs_PDL47_DEGs_DoRothEA.csv")

#------------------------------------------------------------------------------------------------------------
# TF enrichment using TPM and DEGs
# Read object generated above
TPM <- "TPM_normalized_data_cutoff_mean_5_proteincoding.csv"
PDL_24_vs_PDL36 <- "PDL24_vs_PDL36_DEGs_DoRothEA.csv"
PDL_24_vs_PDL47 <- "PDL24_vs_PDL47_DEGs_DoRothEA.csv"
PDL_36_vs_PDL47 <- "PDL36_vs_PDL47_DEGs_DoRothEA.csv"

# Import
PDL_24_vs_PDL36_data <- read_csv(PDL_24_vs_PDL36) 
PDL_24_vs_PDL47_data <- read_csv(PDL_24_vs_PDL47) 
PDL_36_vs_PDL47_data <- read_csv(PDL_36_vs_PDL47) 
TPM<- read.csv(TPM,header=T,stringsAsFactors = F) 

#Change sample same for TPM data
colnames(TPM) <- c("ENSEMBL","external_gene_name","gene_biotype",
                   "PDL24_rep1","PDL24_rep2","PDL24_rep3",
                   "PDL36_rep1","PDL36_rep2","PDL36_rep3",
                   "PDL47_rep1","PDL47_rep2","PDL47_rep3") 
colnames(TPM)
head(TPM)
#Count genes TPM>5
length(TPM$ENSEMBL)
#library("biomaRt")
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Show list of attributes 
listAttributes(mart)

#change names of genes (add names of genes)
gene_list_PDL_24_vs_PDL36 <-getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                                  filters = 'ensembl_gene_id', 
                                  values = PDL_24_vs_PDL36_data,
                                  mart = mart)
gene_list_PDL_24_vs_PDL47 <-getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                                  filters = 'ensembl_gene_id', 
                                  values = PDL_24_vs_PDL47_data,
                                  mart = mart)
gene_list_PDL_36_vs_PDL47 <-getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                                  filters = 'ensembl_gene_id', 
                                  values = PDL_36_vs_PDL47_data,
                                  mart = mart)
#Put SYMBOL
PDL_24_vs_PDL36_data <- inner_join(PDL_24_vs_PDL36_data,gene_list_PDL_24_vs_PDL36,by="ensembl_gene_id")
PDL_24_vs_PDL47_data <- inner_join(PDL_24_vs_PDL47_data,gene_list_PDL_24_vs_PDL47,by="ensembl_gene_id")
PDL_36_vs_PDL47_data <- inner_join(PDL_36_vs_PDL47_data,gene_list_PDL_36_vs_PDL47,by="ensembl_gene_id")

#slightly modify the format of the input files to make them suitable for running DoRothEA.
PDL_24_vs_PDL36_Input <- PDL_24_vs_PDL36_data %>% dplyr::select(external_gene_name,log2FoldChange)%>% dplyr::filter(!is.na(log2FoldChange))%>%column_to_rownames(var="external_gene_name")%>%as.matrix() 
PDL_24_vs_PDL47_Input <- PDL_24_vs_PDL47_data %>% dplyr::select(external_gene_name,log2FoldChange)%>% dplyr::filter(!is.na(log2FoldChange))%>%column_to_rownames(var="external_gene_name")%>%as.matrix() 
PDL_36_vs_PDL47_Input <- PDL_36_vs_PDL47_data %>% dplyr::select(external_gene_name,log2FoldChange)%>% dplyr::filter(!is.na(log2FoldChange))%>%column_to_rownames(var="external_gene_name")%>%as.matrix() 

#TPM
new_name <- make.names(TPM$external_gene_name, unique=TRUE)
rownames(TPM) <- new_name
TPM <- TPM[,c(4:ncol(TPM))]
TPM <- as.matrix(TPM)
TPM

# We load Dorothea Regulons
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

#Run
tf_activities_stat_PDL24_PDL36 <- dorothea::run_viper(PDL_24_vs_PDL36_Input, regulons,
                                                      options =  list(minsize = 5, eset.filter = FALSE, 
                                                                      cores = 1, verbose = FALSE, nes = TRUE))
tf_activities_stat_PDL24_PDL47 <- dorothea::run_viper(PDL_24_vs_PDL47_Input, regulons,
                                                      options =  list(minsize = 5, eset.filter = FALSE, 
                                                                      cores = 1, verbose = FALSE, nes = TRUE))
tf_activities_stat_PDL36_PDL47 <- dorothea::run_viper(PDL_36_vs_PDL47_Input, regulons,
                                                      options =  list(minsize = 5, eset.filter = FALSE, 
                                                                      cores = 1, verbose = FALSE, nes = TRUE))
tf_activities_stat_PDL36_PDL47

#Show each top 10
tf_activities_stat_top10_PDL24_PDL36 <- tf_activities_stat_PDL24_PDL36 %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "log2FoldChange") %>%
  dplyr::top_n(10, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

tf_activities_stat_top10_PDL24_PDL47 <- tf_activities_stat_PDL24_PDL47 %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "log2FoldChange") %>%
  dplyr::top_n(10, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

tf_activities_stat_top10_PDL36_PDL47 <- tf_activities_stat_PDL36_PDL47 %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "log2FoldChange") %>%
  dplyr::top_n(10, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

#Integrate
tf_activities_stat <- full_join(tf_activities_stat_top10_PDL24_PDL36,tf_activities_stat_top10_PDL24_PDL47,by="GeneID")
tf_activities_stat <- full_join(tf_activities_stat,tf_activities_stat_top10_PDL36_PDL47,by="GeneID")
tf_activities_stat
#Import TPM value
tf_activities_counts <- 
  dorothea::run_viper(TPM, regulons,
                      options =  list(minsize = 5, eset.filter = FALSE, 
                                      cores = 1, verbose = FALSE, method = c("scale")))
tf_activities_counts

tf_activities_counts_filter <- tf_activities_counts %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::filter(GeneID %in% tf_activities_stat$GeneID) %>%
  column_to_rownames(var = "GeneID") %>%
  as.matrix()

#Take out vector
tf_activities_vector <- as.vector(tf_activities_counts_filter)

#Make heatmap
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

dorotheaBreaks <- c(seq(min(tf_activities_vector), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(tf_activities_vector)/paletteLength, 
                        max(tf_activities_vector), 
                        length.out=floor(paletteLength/2)))
tf_activities_counts_filter
dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=7, fontsize_row = 7, fontsize_col = 7,fontfamily = "sans", 
                          color=myColor, breaks = dorotheaBreaks,cluster_cols = F, 
                          clustering_distance_rows = "correlation",
                          clustering_method = "ward.D2",
                          main = "Transcription Factor activity", angle_col = 45,
                          treeheight_col = 0,  border_color = NA,
                          width =3,height=3,legend=F,
                          filename = "Fig1_E_drothea_TPM_paper.png")
