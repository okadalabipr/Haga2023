# 10_Heatmap_correlation_invivo_invitro.R
# Make Venn figure between public in vivo and in vitro
# Make heat map of KEGG analysis in vitro and vitro
# Make correlation figure of in vitro and vitro

# Libraies I need
library("tidyverse")
library("edgeR")
library("biomaRt")
library("rrcov")
library("DESeq2")
library("pheatmap")
library("clusterProfiler")
library("org.Hs.eg.db")
library("ggpmisc")
library("export")
library("ggrepel")

# TPM (transcripts per million) calculation for in vivo samples
# Set were to save data
setwd("d:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/TPM/") 

# Check were you are now
getwd() 

# clear the decks
rm(list = ls())

# Make object d (read count data for in vivo)
d <- read.table("d:/haga/human_skin/counts/counts.txt", header=T, stringsAsFactors = F)

# d[,-2:-5]
# remove unwanted lines
d <- d[,!colnames(d) %in% c("Chr", "Start", "End", "Strand")]

# change names of samples(initial names are long)
colnames(d) <- c("Geneid","Length", "SRR7093809", "SRR7093810", "SRR7093811", "SRR7093812", "SRR7093813", "SRR7093814", "SRR7093815", "SRR7093816","SRR7093817", "SRR7093818", "SRR7093819", "SRR7093820", "SRR7093821","SRR7093822", "SRR7093823","SRR7093824", "SRR7093825", "SRR7093826", "SRR7093827", "SRR7093828","SRR7093829", "SRR7093830","SRR7093831", "SRR7093832", "SRR7093833", "SRR7093834", "SRR7093835","SRR7093836","SRR7093837", "SRR7093838","SRR7093839", "SRR7093840", "SRR7093841", "SRR7093842", "SRR7093843","SRR7093844","SRR7093845","SRR7093846", "SRR7093847","SRR7093848", "SRR7093849", "SRR7093850", "SRR7093851", "SRR7093852","SRR7093853","SRR7093854","SRR7093855", "SRR7093856","SRR7093857", "SRR7093858", "SRR7093859", "SRR7093860", "SRR7093861","SRR7093862","SRR7093863","SRR7093864", "SRR7093865","SRR7093866", "SRR7093867", "SRR7093868", "SRR7093869", "SRR7093870","SRR7093871","SRR7093872","SRR7093873", "SRR7093874","SRR7093875", "SRR7093876", "SRR7093877", "SRR7093878", "SRR7093879","SRR7093880","SRR7093881","SRR7093882", "SRR7093883","SRR7093884", "SRR7093885", "SRR7093886", "SRR7093887", "SRR7093888","SRR7093889","SRR7093890","SRR7093891", "SRR7093892","SRR7093893", "SRR7093894", "SRR7093895", "SRR7093896", "SRR7093897","SRR7093898","SRR7093899","SRR7093900", "SRR7093901","SRR7093902", "SRR7093903", "SRR7093904", "SRR7093905", "SRR7093906","SRR7093907", "SRR7093908","SRR7093909","SRR7093910","SRR7093911", "SRR7093912","SRR7093913", "SRR7093914", "SRR7093915", "SRR7093916", "SRR7093917","SRR7093918","SRR7093919", "SRR7093920","SRR7093921", "SRR7093922", "SRR7093923", "SRR7093924", "SRR7093925", "SRR7093926","SRR7093927","SRR7093928", "SRR7093929","SRR7093930", "SRR7093931", "SRR7093932", "SRR7093933", "SRR7093934","SRR7093935","SRR7093936","SRR7093937", "SRR7093938","SRR7093939", "SRR7093940", "SRR7093941", "SRR7093942", "SRR7093943","SRR7093944","SRR7093945","SRR7093946", "SRR7093947","SRR7093948", "SRR7093949", "SRR7093950", "SRR7093951")
d_2 <- d[,c(1,2,
            81,82,83,91,
            93,5,124,6,7,122,8,9,10,11,12,
            13,14,15,16,102,103,
            19,26,24,17,27,18,28,29,22,30,31,20,25,32,
            21,23,104,106,
            33,34,35,36,37,38,39,40,41,42,133,43,44,135,
            134,109,
            131,45,46,47,48,49,50,51,52,53,54,
            56,57,118,58,59,60,115,120,62,63,114,117,
            64,116,65,66,119,113,67)]

colnames(d_2)<- c("ensembl_gene_id", "Length", "SRR7093887_10s","SRR7093888_10s","SRR7093889_10s","SRR7093897_10s","SRR7093899_20s","SRR7093811_20s","SRR7093930_20s","SRR7093812_20s","SRR7093813_20s","SRR7093928_20s","SRR7093814_20s",
                  "SRR7093815_20s","SRR7093816_20s","SRR7093817_20s","SRR7093818_20s","SRR7093819_30s","SRR7093820_30s","SRR7093821_30s","SRR7093822_30s","SRR7093908_30s","SRR7093909_30s","SRR7093825_40s",
                  "SRR7093832_40s","SRR7093830_40s","SRR7093823_40s","SRR7093833_40s","SRR7093824_40s","SRR7093834_40s","SRR7093835_40s","SRR7093828_40s","SRR7093836_40s","SRR7093837_40s","SRR7093826_40s",
                  "SRR7093831_40s","SRR7093838_40s","SRR7093827_50s","SRR7093829_50s","SRR7093910_50s","SRR7093912_50s","SRR7093839_60s","SRR7093840_60s","SRR7093841_60s","SRR7093842_60s","SRR7093843_60s",
                  "SRR7093844_60s","SRR7093845_60s","SRR7093846_60s","SRR7093847_60s","SRR7093848_60s","SRR7093939_60s","SRR7093849_60s","SRR7093850_60s","SRR7093941_60s","SRR7093940_70s","SRR7093915_70s",
                  "SRR7093937_early_80s","SRR7093851_early_80s","SRR7093852_early_80s","SRR7093853_early_80s","SRR7093854_early_80s","SRR7093855_early_80s","SRR7093856_early_80s","SRR7093857_early_80s","SRR7093858_early_80s",
                  "SRR7093859_early_80s","SRR7093860_early_80s","SRR7093862_late_80s","SRR7093863_late_80s","SRR7093924_late_80s","SRR7093864_late_80s","SRR7093865_late_80s","SRR7093866_late_80s","SRR7093921_late_80s",
                  "SRR7093926_late_80s","SRR7093868_late_80s","SRR7093869_late_80s","SRR7093920_late_80s","SRR7093923_late_80s","SRR7093870_90s","SRR7093922_90s","SRR7093871_90s","SRR7093872_90s","SRR7093925_90s",
                  "SRR7093919_90s","SRR7093873_90s")

# Divide the sample read count by the length of the gene multiplied by 1000
d_tpm <- (d_2[,3:ncol(d_2)]/d_2$Length) *1000

# combind Geneid from d to d_tpm
d_tpm <- cbind(ensembl_gene_id = d_2$ensembl_gene_id, d_tpm)

# Get the normalization factor, if group doesn't exist, will be in same group
normfactor<- DGEList(counts = d_tpm[,2:ncol(d_tpm)], group = colnames(d_tpm[,2:ncol(d_tpm)]))
normfactor <- calcNormFactors(normfactor, method="RLE")

# - $counts 
normfactor_samples <- normfactor$samples

# multiply normalization factor with the library size
normfactor_samples$normlib <- normfactor_samples$lib.size*normfactor_samples$norm.factors

# Loop to get the final TPM counts table
for(i in 1:(dim(d_tpm)[2]-1)){
  
  d_tpm[,i+1] <- (d_tpm[,i+1]/normfactor_samples$normlib[i])*1000000
  
}

#Library read
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

#cut off version information
d_tpm[,1] <- str_replace(d_tpm[,1],
                         pattern = ".[0-9]+$",
                         replacement = "")
#chage names of genes (add names of genes)
gene_list <-getBM(attributes = c("ensembl_gene_id","external_gene_name","gene_biotype"),
                  filters = 'ensembl_gene_id', 
                  values = d_tpm,
                  mart = mart)

# Combing gene_list & d_tpm
Final_gene_list <- dplyr::inner_join(gene_list, d_tpm ,by="ensembl_gene_id")
sum(table(Final_gene_list$ensembl_gene_id))
norm_count_mean <- mutate (Final_gene_list, mean = rowMeans(Final_gene_list[,4:ncol(Final_gene_list)]))
norm_count_mean$sig <- if_else(norm_count_mean$mean>5, "TRUE", "FALSE")

# True:11361
# False:49249
# sum:60610
# Coding genes:19937
table(norm_count_mean$sig)
sum(table(norm_count_mean$sig))
table(norm_count_mean$gene_biotype)

# Filter True
norm_count_mean_TRUE <- norm_count_mean %>% dplyr::filter(sig == TRUE)
table(norm_count_mean_TRUE$gene_biotype)
norm_count <- dplyr::select(norm_count_mean_TRUE, -sig & -mean)
sum(table(norm_count$external_gene_name))

# column to row names
norm_count <- column_to_rownames(norm_count, var = "ensembl_gene_id")

# save csv
write.csv(norm_count, file= "20200601_TPM_normalized_data_cutoff_mean_5.csv")

#------------------------------------------------------------------------------------------------------------
#Perform robust principal component analysis
# Set were to save data
setwd("d:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/PCA/") 

# Check were you are now
getwd() 

# clear the decks
rm(list = ls())

# object
counts <- read.csv("D:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/TPM/20200601_TPM_normalized_data_cutoff_mean_5.csv",header = TRUE, stringsAsFactors = FALSE,row.names = 1)
counts <- dplyr::select(counts,- "external_gene_name" & -"gene_biotype")
group <- data.frame(sample = factor(c("SRR7093887_10s",       "SRR7093888_10s",       "SRR7093889_10s",       "SRR7093897_10s",       "SRR7093899_20s",      
                                      "SRR7093811_20s",       "SRR7093930_20s",       "SRR7093812_20s",       "SRR7093813_20s",       "SRR7093928_20s",      
                                      "SRR7093814_20s",       "SRR7093815_20s",       "SRR7093816_20s",       "SRR7093817_20s",       "SRR7093818_20s",      
                                      "SRR7093819_30s",       "SRR7093820_30s",       "SRR7093821_30s",       "SRR7093822_30s",       "SRR7093908_30s",      
                                      "SRR7093909_30s",       "SRR7093825_40s",       "SRR7093832_40s",       "SRR7093830_40s",       "SRR7093823_40s",      
                                      "SRR7093833_40s",       "SRR7093824_40s",       "SRR7093834_40s",       "SRR7093835_40s",       "SRR7093828_40s",      
                                      "SRR7093836_40s",       "SRR7093837_40s",       "SRR7093826_40s",       "SRR7093831_40s",       "SRR7093838_40s",      
                                      "SRR7093827_50s",       "SRR7093829_50s",       "SRR7093910_50s",       "SRR7093912_50s",       "SRR7093839_60s",      
                                      "SRR7093840_60s",       "SRR7093841_60s",       "SRR7093842_60s",       "SRR7093843_60s",       "SRR7093844_60s",      
                                      "SRR7093845_60s",       "SRR7093846_60s",       "SRR7093847_60s",       "SRR7093848_60s",       "SRR7093939_60s",      
                                      "SRR7093849_60s",       "SRR7093850_60s",       "SRR7093941_60s",       "SRR7093940_70s",       "SRR7093915_70s",      
                                      "SRR7093937_early_80s", "SRR7093851_early_80s", "SRR7093852_early_80s", "SRR7093853_early_80s", "SRR7093854_early_80s",
                                      "SRR7093855_early_80s", "SRR7093856_early_80s", "SRR7093857_early_80s", "SRR7093858_early_80s", "SRR7093859_early_80s",
                                      "SRR7093860_early_80s", "SRR7093862_late_80s",  "SRR7093863_late_80s",  "SRR7093924_late_80s" , "SRR7093864_late_80s", 
                                      "SRR7093865_late_80s" , "SRR7093866_late_80s",  "SRR7093921_late_80s",  "SRR7093926_late_80s" , "SRR7093868_late_80s", 
                                      "SRR7093869_late_80s" , "SRR7093920_late_80s",  "SRR7093923_late_80s",  "SRR7093870_90s",       "SRR7093922_90s",      
                                      "SRR7093871_90s",       "SRR7093872_90s",       "SRR7093925_90s",       "SRR7093919_90s",       "SRR7093873_90s"    )),
                    group_name = factor(c("10s","10s","10s","10s",
                                          "20s","20s","20s","20s","20s","20s","20s","20s","20s","20s","20s",
                                          "30s","30s","30s","30s","30s","30s",
                                          "40s","40s","40s","40s","40s","40s","40s","40s","40s","40s","40s","40s","40s","40s",
                                          "50s","50s","50s","50s",
                                          "60s","60s","60s","60s","60s","60s","60s","60s","60s","60s","60s","60s","60s","60s",
                                          "70s","70s",
                                          "early_80s","early_80s","early_80s","early_80s","early_80s","early_80s","early_80s","early_80s","early_80s","early_80s",
                                          "early_80s","late_80s","late_80s","late_80s","late_80s","late_80s","late_80s","late_80s","late_80s","late_80s",
                                          "late_80s","late_80s","late_80s",
                                          "90s","90s","90s","90s","90s","90s","90s")))

counts_t <- t(counts)%>%as.data.frame()                                         
counts_t <- rownames_to_column(counts_t, var ="sample")
counts_group <- inner_join(group, counts_t, by ="sample")

# 10s-20s
counts_group_10_20s <- dplyr::filter(counts_group,counts_group$group_name =="10s"|counts_group$group_name =="20s")
counts_group_10_20s <- column_to_rownames(counts_group_10_20s, var = "sample")

# 30s-40s
counts_group_30s_40s <- dplyr::filter(counts_group,counts_group$group_name =="30s"|counts_group$group_name =="40s")
counts_group_30s_40s <- column_to_rownames(counts_group_30s_40s, var = "sample")

# 60s
counts_group_50s_60s_70s <- dplyr::filter(counts_group,counts_group$group_name =="50s"|counts_group$group_name =="60s"|counts_group$group_name =="70s")
counts_group_50s_60s_70s <- column_to_rownames(counts_group_50s_60s_70s, var = "sample")

# 60s,70s
counts_group_60s_70s <- dplyr::filter(counts_group,counts_group$group_name =="60s"|counts_group$group_name =="70s")
counts_group_60s_70s <- column_to_rownames(counts_group_60s_70s, var = "sample")

# early_80s
counts_group_early_80s <- dplyr::filter(counts_group,counts_group$group_name =="early_80s")
counts_group_early_80s <- column_to_rownames(counts_group_early_80s, var = "sample")

# late_80s
counts_group_late_80s <- dplyr::filter(counts_group,counts_group$group_name =="late_80s")
counts_group_late_80s <- column_to_rownames(counts_group_late_80s, var = "sample")

# 90s
counts_group_90s <- dplyr::filter(counts_group,counts_group$group_name =="90s")
counts_group_90s <- column_to_rownames(counts_group_90s, var = "sample")

#10s-20s
ROBPCA_10s_20s <- PcaHubert(counts_group_10_20s[,2:ncol(counts_group_10_20s)],k=2, scale = FALSE) 
plot(ROBPCA_10s_20s, cex =3) 
outliers_10s_20s <- which(ROBPCA_10s_20s@flag=="TRUE")
write.csv(outliers_10s_20s, file = "20200601_outliers_omit_10s_20s.csv")
#30s-40s
ROBPCA_30s_40s <- PcaHubert(counts_group_30s_40s[,2:ncol(counts_group_30s_40s)],k=2, scale = FALSE) 
plot(ROBPCA_30s_40s, cex =3)
outliers_30s_40s <- which(ROBPCA_30s_40s@flag=="TRUE")
write.csv(outliers_30s_40s, file = "20200601_outliers_omit_30s_40s.csv")

ROBPCA_50s_60s_70s <- PcaHubert(counts_group_50s_60s_70s[,2:ncol(counts_group_50s_60s_70s)],k=2, scale = FALSE) 
plot(ROBPCA_50s_60s_70s, cex =5)
outliers_50s_60s_70s <- which(ROBPCA_50s_60s_70s@flag=="TRUE")
write.csv(outliers_50s_60s_70s, file = "20200601_outliers_omit_50s_60s_70s.csv")

ROBPCA_early_80s <- PcaHubert(counts_group_early_80s[,2:ncol(counts_group_early_80s)],k=2, scale = FALSE) 
plot(ROBPCA_early_80s, cex =5)
outliers_early_80s <- which(ROBPCA_early_80s@flag=="TRUE")
write.csv(outliers_early_80s, file = "20200601_outliers_omit_early80s.csv")

ROBPCA_late_80s <- PcaHubert(counts_group_late_80s[,2:ncol(counts_group_late_80s)],k=2, scale = FALSE) 
plot(ROBPCA_late_80s, cex =5)
outliers_late_80s <- which(ROBPCA_late_80s@flag=="TRUE")
write.csv(outliers_late_80s, file = "20200601_outliers_omit_late80s.csv")

ROBPCA_90s <- PcaHubert(counts_group_90s[,2:ncol(counts_group_90s)],k=2, scale = FALSE) 
plot(ROBPCA_90s, cex =5)
outliers_90s <- which(ROBPCA_90s@flag=="TRUE")
write.csv(outliers_90s, file = "20200601_outliers_omit_90s.csv")

# save outliers
a <- read.csv("20200601_outliers_omit_10s_20s.csv")
b <- read.csv("20200601_outliers_omit_30s_40s.csv")
c <- read.csv("20200601_outliers_omit_50s_60s_70s.csv")
d <- read.csv("20200601_outliers_omit_early80s.csv")
e <- read.csv("20200601_outliers_omit_late80s.csv")
g <- read.csv("20200601_outliers_omit_90s.csv")

outliers_omit <- full_join(a,b,by ="X")
outliers_omit <- full_join(outliers_omit,c,by ="X")
outliers_omit <- full_join(outliers_omit,d,by ="X")
outliers_omit <- full_join(outliers_omit,e,by ="X")
outliers_omit <- full_join(outliers_omit,g,by ="X")

write.csv(outliers, file = "20200601_outliers_omit.csv")

# outliers_omitted TPM
outliers_omit <- dplyr::select(outliers_omit, X)
outliers_omit <- dplyr::rename(outliers_omit, "sample"="X")
counts_outliers <- inner_join(outliers_omit, counts_t, by ="sample") %>% column_to_rownames(var = "sample")
counts_outliers <- t(counts_outliers) %>% as.data.frame()

# save
write.csv(counts_outliers, file = "20200602_TPM_outlier_omitted.csv")

#------------------------------------------------------------------------------------------------------------
# Set were to save data
setwd("d:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/sample_clustering/") 

# Check were you are now
getwd() 

# clear the decks
rm(list = ls())

# read data
norm_count <- read.csv("D:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/PCA/20200602_TPM_outlier_omitted.csv", header = T, stringsAsFactors = F, row.names = 1)

# Select 10s~70s
norm_selected <- norm_count[,c(1:45)]
colnames(norm_selected)
rho <- cor(norm_selected, method = "spearman")

# plot
d <- as.dist(1 - rho)
h <- hclust(d, method = "ward.D2")
plot(h)
normalized_count_data_selected <- rownames_to_column(norm_selected)

# Change row names
names(normalized_count_data_selected)[1] <- "ensembl_gene_id"

# cluster
result <- cutree(h, k=3)
result <- as.data.frame(result)

#Import age info.
age <- read.csv("d:/haga/human_skin/sample/human_skin_sample_data_R.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
age <- dplyr::select(age, Age) 
age
# chage names
rownames(age)<- c("SRR7093887_10s","SRR7093888_10s","SRR7093889_10s","SRR7093897_10s","SRR7093899_20s","SRR7093811_20s","SRR7093930_20s","SRR7093812_20s","SRR7093813_20s","SRR7093928_20s","SRR7093814_20s",
                  "SRR7093815_20s","SRR7093816_20s","SRR7093817_20s","SRR7093818_20s","SRR7093819_30s","SRR7093820_30s","SRR7093821_30s","SRR7093822_30s","SRR7093908_30s","SRR7093909_30s","SRR7093825_40s",
                  "SRR7093832_40s","SRR7093830_40s","SRR7093823_40s","SRR7093833_40s","SRR7093824_40s","SRR7093834_40s","SRR7093835_40s","SRR7093828_40s","SRR7093836_40s","SRR7093837_40s","SRR7093826_40s",
                  "SRR7093831_40s","SRR7093838_40s","SRR7093827_50s","SRR7093829_50s","SRR7093910_50s","SRR7093912_50s","SRR7093839_60s","SRR7093840_60s","SRR7093841_60s","SRR7093842_60s","SRR7093843_60s",
                  "SRR7093844_60s","SRR7093845_60s","SRR7093846_60s","SRR7093847_60s","SRR7093848_60s","SRR7093939_60s","SRR7093849_60s","SRR7093850_60s","SRR7093941_60s","SRR7093940_70s","SRR7093915_70s",
                  "SRR7093937_early_80s","SRR7093851_early_80s","SRR7093852_early_80s","SRR7093853_early_80s","SRR7093854_early_80s","SRR7093855_early_80s","SRR7093856_early_80s","SRR7093857_early_80s","SRR7093858_early_80s",
                  "SRR7093859_early_80s","SRR7093860_early_80s","SRR7093862_late_80s","SRR7093863_late_80s","SRR7093924_late_80s","SRR7093864_late_80s","SRR7093865_late_80s","SRR7093866_late_80s","SRR7093921_late_80s",
                  "SRR7093926_late_80s","SRR7093868_late_80s","SRR7093869_late_80s","SRR7093920_late_80s","SRR7093923_late_80s","SRR7093870_90s","SRR7093922_90s","SRR7093871_90s","SRR7093872_90s","SRR7093925_90s",
                  "SRR7093919_90s","SRR7093873_90s")

# inner_join
age <- rownames_to_column(age)
result <- rownames_to_column(result)
cluster_age <- inner_join(result,age,by ="rowname")

# cluster info
# cluster 1  2  3 
#        22 13 10 
table(cluster_age$result)
cluster_age

#cluster age
age_cluster <- group_by(cluster_age, result)
age_cluster
Cluster1 <- age_cluster %>% dplyr::filter(result=="1")
Cluster2 <- age_cluster %>% dplyr::filter(result=="2")
Cluster3 <- age_cluster %>% dplyr::filter(result=="3")
ave
# p>0.05=Man-Whitney U test
shapiro.test(Cluster1$Age)
shapiro.test(Cluster2$Age)
shapiro.test(Cluster3$Age)

#p value calculation
wilcox.test(Cluster1$Age,Cluster2$Age)
wilcox.test(Cluster2$Age,Cluster3$Age)
wilcox.test(Cluster1$Age,Cluster3$Age)

ave <- age_cluster %>% summarise(mean_age =mean(Age), median_age = median(Age))
ave
# box plot
g <- ggplot(cluster_age, aes(x = result, y = Age, group=result,fill=result, colour=result))
g <- g + geom_boxplot(size=1,outlier.colour = NA,fill="gray75")
g <- g + geom_point(position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.5, jitter.height = 0),alpha = 0.6 , size = 2,color="black")
g <- g + theme_bw()
g <- g+ theme(axis.title = element_text(size = 8), 
              axis.text = element_text(size = 8), 
              axis.text.x = element_text(size = 8), 
              axis.text.y = element_text(size = 8), 
              legend.text = element_text(size = 8), 
              legend.title = element_text(size = 8)) + theme(plot.title = element_text(size = 8)) + xlab("cluster") + ylab("Age distribution in clusters")
g <- g + scale_y_continuous(expand = c(0, 0), limits = c(0,90))
g <- g + theme(legend.position = "none")
g <- g + geom_text(size=8,x = 2, y = 85, label = "*",color="black") 
g <- g + geom_segment(x = 1, xend = 1, y = 75, yend = 80,color="black") 
g <- g + geom_segment(x = 1, xend = 3, y = 80, yend = 80,color="black") 
g <- g + geom_segment(x = 3, xend = 3, y = 75, yend = 80,color="black")
g
#------------------------------------------------------------------------------------------------------------
#Obtain DEGs for in vivo clusters
# Set were to save data
setwd("d:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/DESeq2/") 

# Check were you are now
getwd() 

# clear the decks
rm(list = ls())

# Make object (stringsAsFactors = F means ignore letters)
counts_deseq <- read.table("D:/haga/human_skin/counts/counts.txt", header=T, stringsAsFactors = F)

# remove unwanted lines
counts_deseq <- counts_deseq[,!colnames(counts_deseq) %in% c("Chr", "Start", "End", "Strand")]

# change names of samples(initial names are long)
colnames(counts_deseq) <- c("Geneid", "Lengh", "SRR7093809", "SRR7093810", "SRR7093811", "SRR7093812", "SRR7093813", "SRR7093814", "SRR7093815", "SRR7093816","SRR7093817", "SRR7093818", "SRR7093819", "SRR7093820", "SRR7093821","SRR7093822", "SRR7093823","SRR7093824", "SRR7093825", "SRR7093826", "SRR7093827", "SRR7093828","SRR7093829", "SRR7093830","SRR7093831", "SRR7093832", "SRR7093833", "SRR7093834", "SRR7093835","SRR7093836","SRR7093837", "SRR7093838","SRR7093839", "SRR7093840", "SRR7093841", "SRR7093842", "SRR7093843","SRR7093844","SRR7093845","SRR7093846", "SRR7093847","SRR7093848", "SRR7093849", "SRR7093850", "SRR7093851", "SRR7093852","SRR7093853","SRR7093854","SRR7093855", "SRR7093856","SRR7093857", "SRR7093858", "SRR7093859", "SRR7093860", "SRR7093861","SRR7093862","SRR7093863","SRR7093864", "SRR7093865","SRR7093866", "SRR7093867", "SRR7093868", "SRR7093869", "SRR7093870","SRR7093871","SRR7093872","SRR7093873", "SRR7093874","SRR7093875", "SRR7093876", "SRR7093877", "SRR7093878", "SRR7093879","SRR7093880","SRR7093881","SRR7093882", "SRR7093883","SRR7093884", "SRR7093885", "SRR7093886", "SRR7093887", "SRR7093888","SRR7093889","SRR7093890","SRR7093891", "SRR7093892","SRR7093893", "SRR7093894", "SRR7093895", "SRR7093896", "SRR7093897","SRR7093898","SRR7093899","SRR7093900", "SRR7093901","SRR7093902", "SRR7093903", "SRR7093904", "SRR7093905", "SRR7093906","SRR7093907", "SRR7093908","SRR7093909","SRR7093910","SRR7093911", "SRR7093912","SRR7093913", "SRR7093914", "SRR7093915", "SRR7093916", "SRR7093917","SRR7093918","SRR7093919", "SRR7093920","SRR7093921", "SRR7093922", "SRR7093923", "SRR7093924", "SRR7093925", "SRR7093926","SRR7093927","SRR7093928", "SRR7093929","SRR7093930", "SRR7093931", "SRR7093932", "SRR7093933", "SRR7093934","SRR7093935","SRR7093936","SRR7093937", "SRR7093938","SRR7093939", "SRR7093940", "SRR7093941", "SRR7093942", "SRR7093943","SRR7093944","SRR7093945","SRR7093946", "SRR7093947","SRR7093948", "SRR7093949", "SRR7093950", "SRR7093951")
counts_deseq <- counts_deseq[,c(1,
                                81,82,83,91,
                                93,5,124,6,7,122,8,9,10,11,12,
                                13,14,15,16,102,103,
                                19,26,24,17,27,18,28,29,22,30,31,20,25,32,
                                21,23,104,106,
                                33,34,35,36,37,38,39,40,41,42,133,43,44,135,
                                134,109,
                                131,45,46,47,48,49,50,51,52,53,54,
                                56,57,118,58,59,60,115,120,62,63,114,117,
                                64,116,65,66,119,113,67)]

# Change names
colnames(counts_deseq)<- c("ensembl_gene_id","SRR7093887_10s","SRR7093888_10s","SRR7093889_10s","SRR7093897_10s","SRR7093899_20s","SRR7093811_20s","SRR7093930_20s","SRR7093812_20s","SRR7093813_20s","SRR7093928_20s","SRR7093814_20s",
                           "SRR7093815_20s","SRR7093816_20s","SRR7093817_20s","SRR7093818_20s","SRR7093819_30s","SRR7093820_30s","SRR7093821_30s","SRR7093822_30s","SRR7093908_30s","SRR7093909_30s","SRR7093825_40s",
                           "SRR7093832_40s","SRR7093830_40s","SRR7093823_40s","SRR7093833_40s","SRR7093824_40s","SRR7093834_40s","SRR7093835_40s","SRR7093828_40s","SRR7093836_40s","SRR7093837_40s","SRR7093826_40s",
                           "SRR7093831_40s","SRR7093838_40s","SRR7093827_50s","SRR7093829_50s","SRR7093910_50s","SRR7093912_50s","SRR7093839_60s","SRR7093840_60s","SRR7093841_60s","SRR7093842_60s","SRR7093843_60s",
                           "SRR7093844_60s","SRR7093845_60s","SRR7093846_60s","SRR7093847_60s","SRR7093848_60s","SRR7093939_60s","SRR7093849_60s","SRR7093850_60s","SRR7093941_60s","SRR7093940_70s","SRR7093915_70s",
                           "SRR7093937_early_80s","SRR7093851_early_80s","SRR7093852_early_80s","SRR7093853_early_80s","SRR7093854_early_80s","SRR7093855_early_80s","SRR7093856_early_80s","SRR7093857_early_80s","SRR7093858_early_80s",
                           "SRR7093859_early_80s","SRR7093860_early_80s","SRR7093862_late_80s","SRR7093863_late_80s","SRR7093924_late_80s","SRR7093864_late_80s","SRR7093865_late_80s","SRR7093866_late_80s","SRR7093921_late_80s",
                           "SRR7093926_late_80s","SRR7093868_late_80s","SRR7093869_late_80s","SRR7093920_late_80s","SRR7093923_late_80s","SRR7093870_90s","SRR7093922_90s","SRR7093871_90s","SRR7093872_90s","SRR7093925_90s",
                           "SRR7093919_90s","SRR7093873_90s")

#cut off version information
counts_deseq[,1] <- str_replace(counts_deseq[,1],
                                pattern = ".[0-9]+$",
                                replacement = "")

# import data
gene_set <- read.csv("D:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/PCA/20200602_TPM_outlier_omitted.csv", header = T, stringsAsFactors = F)
gene_set <- gene_set$X %>% as.data.frame()
gene_set <- dplyr::rename(gene_set, ensembl_gene_id = .)

# gene_set(TPM>5)
#11361
sum(table(gene_set))
counts_deseq_gene <- dplyr::inner_join(gene_set, counts_deseq, by = "ensembl_gene_id")
Final_gene_list_t <- column_to_rownames(counts_deseq_gene, var ="ensembl_gene_id") %>% t() %>%as.data.frame()
Final_gene_list_t <- rownames_to_column(Final_gene_list_t, var = "X")

# import cluster result
Cluster_result <- read.csv("D:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/sample_clustering/20200602_Cluster_dendrogram_k=3_10s_70s.csv", header = TRUE, stringsAsFactors = FALSE) 

# cluster result 
Final_gene_list_t_1 <- dplyr::inner_join(Cluster_result, Final_gene_list_t,by = "X")
glimpse(Final_gene_list_t_1)

cluster_1_2 <- dplyr::filter(Final_gene_list_t_1, result !=3) 
cluster_1_3 <- dplyr::filter(Final_gene_list_t_1, result !=2) 
cluster_2_3 <- dplyr::filter(Final_gene_list_t_1, result !=1) 

# arrange
cluster_1_2 <- arrange(cluster_1_2, result)
cluster_1_3 <- arrange(cluster_1_3, result)
cluster_2_3 <- arrange(cluster_2_3, result)

#column_to_rownames
cluster_1_2 <- column_to_rownames(cluster_1_2, var = "X")
cluster_1_3 <- column_to_rownames(cluster_1_3, var = "X")
cluster_2_3 <- column_to_rownames(cluster_2_3, var = "X")

cluster_1_2 <- dplyr::select(cluster_1_2, -result) 
cluster_1_3 <- dplyr::select(cluster_1_3, -result) 
cluster_2_3 <- dplyr::select(cluster_2_3, -result) 

cluster_1_2_t <- t(cluster_1_2) %>% as.data.frame()
cluster_1_3_t <- t(cluster_1_3) %>% as.data.frame() 
cluster_2_3_t <- t(cluster_2_3) %>% as.data.frame()


# cluster1_2:35
# cluster1_3:32
# cluster2_3:23

# save.csv
write.csv(cluster_1_2_t, file ="20200602_cluster_1_2.csv")
write.csv(cluster_1_3_t, file ="20200602_cluster_1_3.csv")
write.csv(cluster_2_3_t, file ="20200602_cluster_2_3.csv")

#------------------------------------------------------------------------------------------------------------
# clear desk
rm(list=ls())

# import data above
A <- read.csv("D:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/DESeq2/20200602_cluster_1_2.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
B <- read.csv("D:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/DESeq2/20200602_cluster_1_3.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
C <- read.csv("D:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/DESeq2/20200602_cluster_2_3.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

#Make Tag for each groups
group_A <- data.frame(con = factor(c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1",
                                     "2","2","2","2","2","2","2","2","2","2","2","2","2")))

group_B <- data.frame(con = factor(c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1",
                                     "3","3","3","3","3","3","3","3","3","3")))

group_C <- data.frame(con = factor(c("2","2","2","2","2","2","2","2","2","2","2","2","2",
                                     "3","3","3","3","3","3","3","3","3","3")))

colnames(C)

# Change row name
rownames(group_A) <- c("SRR7093888_10s", "SRR7093889_10s", "SRR7093897_10s", "SRR7093899_20s", "SRR7093811_20s", "SRR7093930_20s",
                       "SRR7093812_20s", "SRR7093928_20s", "SRR7093814_20s", "SRR7093819_30s", "SRR7093908_30s", "SRR7093830_40s",
                       "SRR7093833_40s", "SRR7093834_40s", "SRR7093826_40s", "SRR7093838_40s", "SRR7093829_50s", "SRR7093840_60s",
                       "SRR7093842_60s", "SRR7093846_60s", "SRR7093850_60s", "SRR7093941_60s", "SRR7093813_20s", "SRR7093815_20s",
                       "SRR7093816_20s", "SRR7093817_20s", "SRR7093820_30s", "SRR7093822_30s", "SRR7093835_40s", "SRR7093837_40s",
                       "SRR7093839_60s", "SRR7093841_60s", "SRR7093844_60s", "SRR7093845_60s", "SRR7093848_60s")

rownames(group_B) <- c( "SRR7093888_10s", "SRR7093889_10s", "SRR7093897_10s", "SRR7093899_20s", "SRR7093811_20s", "SRR7093930_20s",
                        "SRR7093812_20s", "SRR7093928_20s", "SRR7093814_20s", "SRR7093819_30s", "SRR7093908_30s", "SRR7093830_40s",
                        "SRR7093833_40s", "SRR7093834_40s", "SRR7093826_40s", "SRR7093838_40s", "SRR7093829_50s", "SRR7093840_60s",
                        "SRR7093842_60s", "SRR7093846_60s", "SRR7093850_60s", "SRR7093941_60s", "SRR7093825_40s", "SRR7093823_40s",
                        "SRR7093824_40s", "SRR7093828_40s", "SRR7093831_40s", "SRR7093827_50s", "SRR7093843_60s", "SRR7093849_60s",
                        "SRR7093940_70s", "SRR7093915_70s")

rownames(group_C) <- c("SRR7093813_20s", "SRR7093815_20s", "SRR7093816_20s", "SRR7093817_20s", "SRR7093820_30s", "SRR7093822_30s",
                       "SRR7093835_40s", "SRR7093837_40s", "SRR7093839_60s", "SRR7093841_60s", "SRR7093844_60s", "SRR7093845_60s",
                       "SRR7093848_60s", "SRR7093825_40s", "SRR7093823_40s", "SRR7093824_40s", "SRR7093828_40s", "SRR7093831_40s",
                       "SRR7093827_50s", "SRR7093843_60s", "SRR7093849_60s", "SRR7093940_70s", "SRR7093915_70s")


# Make_data_set_for_deseq
dds_A <- DESeqDataSetFromMatrix(countData = A, colData = group_A, design = ~ con)
dds_B <- DESeqDataSetFromMatrix(countData = B, colData = group_B, design = ~ con)
dds_C <- DESeqDataSetFromMatrix(countData = C, colData = group_C, design = ~ con)

# Set ref
dds_A$con <- relevel(dds_A$con, ref = "1")
dds_B$con <- relevel(dds_B$con, ref = "1")
dds_C$con <- relevel(dds_C$con, ref = "2")


# check
relevel(dds_A$con, ref = "1")
relevel(dds_B$con, ref = "1")
relevel(dds_C$con, ref = "2")


# Run_deseq
dds_A <- DESeq(dds_A)
dds_B <- DESeq(dds_B)
dds_C <- DESeq(dds_C)

# result
res_A <- results(dds_A)
res_B <- results(dds_B)
res_C <- results(dds_C)

# Make dataset for p value
res_A<- res_A[,c(2, 5)] %>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")
res_B<- res_B[,c(2, 5)] %>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")
res_C<- res_C[,c(2, 5)] %>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")

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

# Caculate adj p value among samples
q <- p.adjust(c(res_A_p,res_B_p,res_C_p), method = "BH") %>% as.data.frame()

# extract each adj p value
q_A <- q %>% dplyr::slice(1:11361) %>% as.data.frame()
q_B <- q %>% dplyr::slice(11362:22722) %>% as.data.frame()
q_C <- q %>% dplyr::slice(22723:nrow(q)) %>% as.data.frame()

#check
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

# Count data of padj<0.05(log2FC<-0.58 or log2FC>0.58, p<0.05)
# result 457
res_A[(res_A$log2FoldChange < -log2(1.5) & res_A$padj < 0.05)  | (res_A$log2FoldChange > log2(1.5) & res_A$padj < 0.05),] %>% nrow()

# Count data of padj<0.05(log2FC<-0.58 or log2FC>0.58, p<0.05)
# result 1545
res_B[(res_B$log2FoldChange < -log2(1.5) & res_B$padj < 0.05)  | (res_B$log2FoldChange > log2(1.5) & res_B$padj < 0.05),] %>% nrow()

# Count data of padj<0.05(log2FC<-0.58 or log2FC>0.58, p<0.05)
# result  1467
res_C[(res_C$log2FoldChange < -log2(1.5) & res_C$padj < 0.05)  | (res_C$log2FoldChange > log2(1.5) & res_C$padj < 0.05),] %>% nrow()

# identify data of p<0.05
res_A$sig <- if_else((res_A$log2FoldChange >log2(1.5) & res_A$padj <0.05) | (res_A$log2FoldChange < -log2(1.5) & res_A$padj <0.05), "TRUE", "FALSE")
res_B$sig <- if_else((res_B$log2FoldChange >log2(1.5) & res_B$padj <0.05) | (res_B$log2FoldChange < -log2(1.5) & res_B$padj <0.05), "TRUE", "FALSE")
res_C$sig <- if_else((res_C$log2FoldChange >log2(1.5) & res_C$padj <0.05) | (res_C$log2FoldChange < -log2(1.5) & res_C$padj <0.05), "TRUE", "FALSE")


#save
#Where: d:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/DESeq2/
write.csv(as.data.frame(res_A), file="20221209_cluster_1_2_deseq_result.csv")
write.csv(as.data.frame(res_B), file="20221209_cluster_1_3_deseq_result.csv")
write.csv(as.data.frame(res_C), file="20221209_cluster_2_3_deseq_result.csv")


#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
# public in vitro data analysis
# TPM (transcripts per million) calculation
# Set were to save data
setwd("k:/haga/skin_vitro/data_analysis_R/TPM/") 

#Project path
Project_path <- c("k:/haga/skin_vitro/data_analysis_R/")
Project_path

# Check were you are now
getwd() 

# clear the decks
rm(list = ls())

# Make object d (stringsAsFactors = F means ignore letters)
d <- read.table("k:/haga/skin_vitro/featureCounts/20201030_output_counts.txt", header=T, stringsAsFactors = F)

# d[,-2:-5]
# remove unwanted lines
d <- d[,!colnames(d) %in% c("Chr", "Start", "End", "Strand")]

# change names of samples(initial names are long)
colnames(d) <- c("Geneid","Length", "SRR1660537", "SRR1660538", "SRR1660539","SRR1660540", "SRR1660541","SRR1660542",
                 "SRR1660543","SRR1660544","SRR1660545","SRR1660546","SRR1660547","SRR1660548",
                 "SRR2751110","SRR2751111","SRR2751112","SRR2751113","SRR2751114","SRR2751115","SRR2751116","SRR2751117","SRR2751118")

colnames(d) <- c("Geneid","Length", "BJ_PD34_1", "BJ_PD34_2", "BJ_PD34_3","BJ_PD72_1", "BJ_PD72_2", "BJ_PD72_3",
                 "HFF_PD16_1","HFF_PD16_2","HFF_PD16_3","HFF_PD74_1","HFF_PD74_2","HFF_PD74_3",
                 "HFF_PD26_1","HFF_PD26_2","HFF_PD26_3","HFF_PD46_1","HFF_PD46_2","HFF_PD46_3","HFF_PD64_1","HFF_PD64_2","HFF_PD64_3")

# Devide with BJ and HFF-1
d_2 <- d[,c(1:11,15:23,12:14)]

# Divide the sample read count by the length of the gene multiplied by 1000
d_tpm <- (d_2[,3:ncol(d_2)]/d_2$Length) *1000

# combind Geneid from d to d_tpm
d_tpm <- cbind(ensembl_gene_id = d_2$Geneid, d_tpm)

# Get the normalization factor, if group doesn't exist, will be in same group
normfactor<- DGEList(counts = d_tpm[,2:ncol(d_tpm)], group = colnames(d_tpm[,2:ncol(d_tpm)]))

normfactor <- calcNormFactors(normfactor, method="RLE")

# - $counts 
normfactor_samples <- normfactor$samples

#write graph
g <- ggplot(normfactor_samples, aes(x =lib.size, y = norm.factors, label=group)) + geom_point(size =4, alpha =0.5) + theme_bw() + xlab("lib.size")+ ylab("norm.factors") + theme(axis.title = element_text(size = 16))
g <- g + scale_y_continuous(limits = c(0, NA))
g <- g + geom_text()
g

# multiply normalization factor with the library size
normfactor_samples$normlib <- normfactor_samples$lib.size*normfactor_samples$norm.factors

# Loop to get the final TPM counts table

for(i in 1:(dim(d_tpm)[2]-1)){
  
  d_tpm[,i+1] <- (d_tpm[,i+1]/normfactor_samples$normlib[i])*1000000
  
}

#Library read
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")


#cut off version information
d_tpm[,1] <- str_replace(d_tpm[,1],
                         pattern = ".[0-9]+$",
                         replacement = "")
#chage names of genes (add names of genes)
gene_list <-getBM(attributes = c("ensembl_gene_id","external_gene_name","gene_biotype"),
                  filters = 'ensembl_gene_id', 
                  values = d_tpm,
                  mart = mart)

# protein_coding; 19912
table(gene_list$gene_biotype)

# Combing gene_list & d_tpm
Final_gene_list <- dplyr::inner_join(gene_list, d_tpm ,by="ensembl_gene_id")
Final_gene_list
sum(table(Final_gene_list$ensembl_gene_id))

norm_count_mean <- mutate (Final_gene_list, mean = rowMeans(Final_gene_list[,4:ncol(Final_gene_list)]))

norm_count_mean$sig <- if_else(norm_count_mean$mean>5, "TRUE", "FALSE")

norm_count_mean
# True:12296
# False:48268
# sum:60564
table(norm_count_mean$sig)
sum(table(norm_count_mean$sig))
table(norm_count_mean$gene_biotype)

# Filter True
norm_count_mean_TRUE <- norm_count_mean %>% dplyr::filter(sig == TRUE)

#Protein conding 10456
table(norm_count_mean_TRUE$gene_biotype)
norm_count <- dplyr::select(norm_count_mean_TRUE, -sig & -mean)
sum(table(norm_count$external_gene_name))

# column to row names
norm_count <- column_to_rownames(norm_count, var = "ensembl_gene_id")

# save csv
write.csv(norm_count, file= "20201030_TPM_normalized_data_cutoff_mean_5.csv")
norm_count

#------------------------------------------------------------------------------------------------------------
# Set were to save data
setwd("k:/haga/skin_vitro/data_analysis_R/DESeq2/") 

# Check were you are now
getwd() 

# clear the decks
rm(list = ls())

# Make object (stringsAsFactors = F means ignore letters)
counts_deseq <- read.table("k:/haga/skin_vitro/featureCounts/20201030_output_counts.txt", header=T, stringsAsFactors = F)

# remove unwanted lines
counts_deseq <- counts_deseq[,!colnames(counts_deseq) %in% c("Chr", "Start", "End", "Strand")]

# change names of samples(initial names are long)
colnames(counts_deseq) <- c("Geneid","Length", "BJ_PD34_1", "BJ_PD34_2", "BJ_PD34_3","BJ_PD72_1", "BJ_PD72_2", "BJ_PD72_3",
                            "HFF_PD16_1","HFF_PD16_2","HFF_PD16_3","HFF_PD74_1","HFF_PD74_2","HFF_PD74_3",
                            "HFF_PD26_1","HFF_PD26_2","HFF_PD26_3","HFF_PD46_1","HFF_PD46_2","HFF_PD46_3","HFF_PD64_1","HFF_PD64_2","HFF_PD64_3")


#cut off version information
counts_deseq[,1] <- str_replace(counts_deseq[,1],
                                pattern = ".[0-9]+$",
                                replacement = "")

counts_deseq <- counts_deseq[,c(1,3:11,15:23,12:14)]
colnames(counts_deseq)

# import data
gene_set <- read.csv("k:/haga/skin_vitro/data_analysis_R/TPM/20201030_TPM_normalized_data_cutoff_mean_5.csv", header = T, stringsAsFactors = F)
gene_set <- gene_set$X %>% as.data.frame()
gene_set <- dplyr::rename(gene_set, Geneid = .)
gene_set

# gene_set(TPM>5)
sum(table(gene_set))
counts_deseq <- dplyr::inner_join(gene_set, counts_deseq, by = "Geneid")
colnames(counts_deseq)

#group
PD16_PD26 <-counts_deseq[,c(1,8:13)]%>% column_to_rownames(var = "Geneid")
PD16_PD46 <-counts_deseq[,c(1,8:10,14:16)]%>% column_to_rownames(var = "Geneid")
PD16_PD64 <-counts_deseq[,c(1,8:10,17:19)]%>% column_to_rownames(var = "Geneid")
PD16_PD74 <-counts_deseq[,c(1,8:10,20:22)]%>% column_to_rownames(var = "Geneid")
PD26_PD46 <-counts_deseq[,c(1,11:13,14:16)]%>% column_to_rownames(var = "Geneid")
PD26_PD64 <-counts_deseq[,c(1,11:13,17:19)]%>% column_to_rownames(var = "Geneid")
PD26_PD74 <-counts_deseq[,c(1,11:13,20:22)]%>% column_to_rownames(var = "Geneid")
PD46_PD64 <-counts_deseq[,c(1,14:16,17:19)]%>% column_to_rownames(var = "Geneid")
PD46_PD74 <-counts_deseq[,c(1,14:16,20:22)]%>% column_to_rownames(var = "Geneid")
PD64_PD74 <-counts_deseq[,c(1,17:19,20:22)]%>% column_to_rownames(var = "Geneid")

#Make Tag for each groups
group_PD16_PD26 <- data.frame(con = factor(c("PD16","PD16","PD16","PD26","PD26","PD26")))
group_PD16_PD46 <- data.frame(con = factor(c("PD16","PD16","PD16","PD46","PD46","PD46")))
group_PD16_PD64 <- data.frame(con = factor(c("PD16","PD16","PD16","PD64","PD64","PD64")))
group_PD16_PD74 <- data.frame(con = factor(c("PD16","PD16","PD16","PD74","PD74","PD74")))
group_PD26_PD46 <- data.frame(con = factor(c("PD26","PD26","PD26","PD46","PD46","PD46")))
group_PD26_PD64 <- data.frame(con = factor(c("PD26","PD26","PD26","PD64","PD64","PD64")))
group_PD26_PD74 <- data.frame(con = factor(c("PD26","PD26","PD26","PD74","PD74","PD74")))
group_PD46_PD64 <- data.frame(con = factor(c("PD46","PD46","PD46","PD64","PD64","PD64")))
group_PD46_PD74 <- data.frame(con = factor(c("PD46","PD46","PD46","PD74","PD74","PD74")))
group_PD64_PD74 <- data.frame(con = factor(c("PD64","PD64","PD64","PD74","PD74","PD74")))

# Change row name
rownames(group_PD16_PD26) <- c("HFF_PD16_1", "HFF_PD16_2", "HFF_PD16_3","HFF_PD26_1","HFF_PD26_2","HFF_PD26_3")
rownames(group_PD16_PD46) <- c("HFF_PD16_1", "HFF_PD16_2", "HFF_PD16_3","HFF_PD46_1","HFF_PD46_2","HFF_PD46_3")
rownames(group_PD16_PD64) <- c("HFF_PD16_1", "HFF_PD16_2", "HFF_PD16_3","HFF_PD64_1","HFF_PD64_2","HFF_PD64_3")
rownames(group_PD16_PD74) <- c("HFF_PD16_1", "HFF_PD16_2", "HFF_PD16_3","HFF_PD74_1","HFF_PD74_2","HFF_PD74_3")

rownames(group_PD26_PD46) <- c("HFF_PD26_1", "HFF_PD26_2", "HFF_PD26_3","HFF_PD46_1","HFF_PD46_2","HFF_PD46_3")
rownames(group_PD26_PD64) <- c("HFF_PD26_1", "HFF_PD26_2", "HFF_PD26_3","HFF_PD64_1","HFF_PD64_2","HFF_PD64_3")
rownames(group_PD26_PD74) <- c("HFF_PD26_1", "HFF_PD26_2", "HFF_PD26_3","HFF_PD74_1","HFF_PD74_2","HFF_PD74_3")

rownames(group_PD46_PD64) <- c("HFF_PD46_1", "HFF_PD46_2", "HFF_PD46_3","HFF_PD64_1","HFF_PD64_2","HFF_PD64_3")
rownames(group_PD46_PD74) <- c("HFF_PD46_1", "HFF_PD46_2", "HFF_PD46_3","HFF_PD74_1","HFF_PD74_2","HFF_PD74_3")

rownames(group_PD64_PD74) <- c("HFF_PD64_1", "HFF_PD64_2", "HFF_PD64_3","HFF_PD74_1","HFF_PD74_2","HFF_PD74_3")

# Make_data_set_for_deseq
dds_2 <- DESeqDataSetFromMatrix(countData = PD16_PD26, colData = group_PD16_PD26, design = ~ con)
dds_3 <- DESeqDataSetFromMatrix(countData = PD16_PD46, colData = group_PD16_PD46, design = ~ con)
dds_4 <- DESeqDataSetFromMatrix(countData = PD16_PD64, colData = group_PD16_PD64, design = ~ con)
dds_5 <- DESeqDataSetFromMatrix(countData = PD16_PD74, colData = group_PD16_PD74, design = ~ con)
dds_6 <- DESeqDataSetFromMatrix(countData = PD26_PD46, colData = group_PD26_PD46, design = ~ con)
dds_7 <- DESeqDataSetFromMatrix(countData = PD26_PD64, colData = group_PD26_PD64, design = ~ con)
dds_8 <- DESeqDataSetFromMatrix(countData = PD26_PD74, colData = group_PD26_PD74, design = ~ con)
dds_9 <- DESeqDataSetFromMatrix(countData = PD46_PD64, colData = group_PD46_PD64, design = ~ con)
dds_10 <- DESeqDataSetFromMatrix(countData = PD46_PD74, colData = group_PD46_PD74, design = ~ con)
dds_11 <- DESeqDataSetFromMatrix(countData = PD64_PD74, colData = group_PD64_PD74, design = ~ con)

# Set ref
dds_2$con <- relevel(dds_2$con, ref = "PD16")
dds_3$con <- relevel(dds_3$con, ref = "PD16")
dds_4$con <- relevel(dds_4$con, ref = "PD16")
dds_5$con <- relevel(dds_5$con, ref = "PD16")
dds_6$con <- relevel(dds_6$con, ref = "PD26")
dds_7$con <- relevel(dds_7$con, ref = "PD26")
dds_8$con <- relevel(dds_8$con, ref = "PD26")
dds_9$con <- relevel(dds_9$con, ref = "PD46")
dds_10$con <- relevel(dds_10$con, ref = "PD46")
dds_11$con <- relevel(dds_11$con, ref = "PD64")

# check
relevel(dds_2$con, ref = "PD16")
relevel(dds_3$con, ref = "PD16")
relevel(dds_4$con, ref = "PD16")
relevel(dds_5$con, ref = "PD16")
relevel(dds_6$con, ref = "PD26")
relevel(dds_7$con, ref = "PD26")
relevel(dds_8$con, ref = "PD26")
relevel(dds_9$con, ref = "PD46")
relevel(dds_10$con, ref = "PD46")
relevel(dds_11$con, ref = "PD64")

# Run_deseq
dds_2 <- DESeq(dds_2)
dds_3 <- DESeq(dds_3)
dds_4 <- DESeq(dds_4)
dds_5 <- DESeq(dds_5)
dds_6 <- DESeq(dds_6)
dds_7 <- DESeq(dds_7)
dds_8 <- DESeq(dds_8)
dds_9 <- DESeq(dds_9)
dds_10 <- DESeq(dds_10)
dds_11 <- DESeq(dds_11)

# result
res_B <- results(dds_2)
res_C <- results(dds_3)
res_D <- results(dds_4)
res_E <- results(dds_5)
res_F <- results(dds_6)
res_G <- results(dds_7)
res_H <- results(dds_8)
res_I <- results(dds_9)
res_J <- results(dds_10)
res_K <- results(dds_11)

# Make dataset for p value
res_B<- res_B[,c(2, 5)]%>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")
res_C<- res_C[,c(2, 5)]%>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")
res_D<- res_D[,c(2, 5)]%>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")
res_E<- res_E[,c(2, 5)]%>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")
res_F<- res_F[,c(2, 5)]%>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")
res_G<- res_G[,c(2, 5)]%>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")
res_H<- res_H[,c(2, 5)]%>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")
res_I<- res_I[,c(2, 5)]%>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")
res_J<- res_J[,c(2, 5)]%>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")
res_K<- res_K[,c(2, 5)]%>% as.data.frame() %>% na.omit()%>% rownames_to_column(var ="ENSEMBL")

#merge files into one
res <- inner_join(res_B,res_C, by = "ENSEMBL")
res <- inner_join(res,res_D, by = "ENSEMBL")
res <- inner_join(res,res_E, by = "ENSEMBL")
res <- inner_join(res,res_F, by = "ENSEMBL")
res <- inner_join(res,res_G, by = "ENSEMBL")
res <- inner_join(res,res_H, by = "ENSEMBL")
res <- inner_join(res,res_I, by = "ENSEMBL")
res <- inner_join(res,res_J, by = "ENSEMBL")
res <- inner_join(res,res_K, by = "ENSEMBL")
colnames(res)
colnames(res) <- c("ENSEMBL",
                   "log2FoldChange_B","pvalue_B",
                   "log2FoldChange_C","pvalue_C",
                   "log2FoldChange_D","pvalue_D",
                   "log2FoldChange_E","pvalue_E",
                   "log2FoldChange_F","pvalue_F",
                   "log2FoldChange_G","pvalue_G",
                   "log2FoldChange_H","pvalue_H",
                   "log2FoldChange_I","pvalue_I",
                   "log2FoldChange_J","pvalue_J",
                   "log2FoldChange_K","pvalue_K")
# adjust p value among samples
res_B_p <- res$pvalue_B
res_C_p <- res$pvalue_C
res_D_p <- res$pvalue_D
res_E_p <- res$pvalue_E
res_F_p <- res$pvalue_F
res_G_p <- res$pvalue_G
res_H_p <- res$pvalue_H
res_I_p <- res$pvalue_I
res_J_p <- res$pvalue_J
res_K_p <- res$pvalue_K
length(res_C_p)

# Caculate  among samples
q <- p.adjust(c(res_B_p,res_C_p,res_D_p,res_E_p,res_F_p,res_G_p,res_H_p,res_I_p,res_J_p,res_K_p), method = "BH") %>% as.data.frame()

# extract each FDR
q_B <- q %>% dplyr::slice(1:12291) %>% as.data.frame()
q_C <- q %>% dplyr::slice(12292:24582) %>% as.data.frame()
q_D <- q %>% dplyr::slice(24583:36873) %>% as.data.frame()
q_E <- q %>% dplyr::slice(36874:49164) %>% as.data.frame()
q_F <- q %>% dplyr::slice(49165:61455) %>% as.data.frame()
q_G <- q %>% dplyr::slice(61456:73746) %>% as.data.frame()
q_H <- q %>% dplyr::slice(73747:86037) %>% as.data.frame()
q_I <- q %>% dplyr::slice(86038:98328) %>% as.data.frame()
q_J <- q %>% dplyr::slice(98329:110619) %>% as.data.frame()
q_K <- q %>% dplyr::slice(110620:nrow(q)) %>% as.data.frame()

#Integrate into one
res_B <- cbind(res$ENSEMBL,res$log2FoldChange_B,q_B) %>% dplyr::rename("padj" = ".")
res_C <- cbind(res$ENSEMBL,res$log2FoldChange_C,q_C) %>% dplyr::rename("padj" = ".")
res_D <- cbind(res$ENSEMBL,res$log2FoldChange_D,q_D) %>% dplyr::rename("padj" = ".")
res_E <- cbind(res$ENSEMBL,res$log2FoldChange_E,q_E) %>% dplyr::rename("padj" = ".")
res_F <- cbind(res$ENSEMBL,res$log2FoldChange_F,q_F) %>% dplyr::rename("padj" = ".")
res_G <- cbind(res$ENSEMBL,res$log2FoldChange_G,q_G) %>% dplyr::rename("padj" = ".")
res_H <- cbind(res$ENSEMBL,res$log2FoldChange_H,q_H) %>% dplyr::rename("padj" = ".")
res_I <- cbind(res$ENSEMBL,res$log2FoldChange_I,q_I) %>% dplyr::rename("padj" = ".")
res_J <- cbind(res$ENSEMBL,res$log2FoldChange_J,q_J) %>% dplyr::rename("padj" = ".")
res_K <- cbind(res$ENSEMBL,res$log2FoldChange_K,q_K) %>% dplyr::rename("padj" = ".")

#rename
colnames(res_B) <- c("ENSEMBL","log2FoldChange","padj")
colnames(res_C) <- c("ENSEMBL","log2FoldChange","padj")
colnames(res_D) <- c("ENSEMBL","log2FoldChange","padj")
colnames(res_E) <- c("ENSEMBL","log2FoldChange","padj")
colnames(res_F) <- c("ENSEMBL","log2FoldChange","padj")
colnames(res_G) <- c("ENSEMBL","log2FoldChange","padj")
colnames(res_H) <- c("ENSEMBL","log2FoldChange","padj")
colnames(res_I) <- c("ENSEMBL","log2FoldChange","padj")
colnames(res_J) <- c("ENSEMBL","log2FoldChange","padj")
colnames(res_K) <- c("ENSEMBL","log2FoldChange","padj")

# identify data of p<0.01
res_B$sig <- if_else((res_B$log2FoldChange >1 & res_B$padj <0.01) | (res_B$log2FoldChange < -1 & res_B$padj <0.01), "TRUE", "FALSE")

res_C$sig <- if_else((res_C$log2FoldChange >1 & res_C$padj <0.01) | (res_C$log2FoldChange < -1 & res_C$padj <0.01), "TRUE", "FALSE")

res_D$sig <- if_else((res_D$log2FoldChange >1 & res_D$padj <0.01) | (res_D$log2FoldChange < -1 & res_D$padj <0.01), "TRUE", "FALSE")

res_E$sig <- if_else((res_E$log2FoldChange >1 & res_E$padj <0.01) | (res_E$log2FoldChange < -1 & res_E$padj <0.01), "TRUE", "FALSE")

res_F$sig <- if_else((res_F$log2FoldChange >1 & res_F$padj <0.01) | (res_F$log2FoldChange < -1 & res_F$padj <0.01), "TRUE", "FALSE")

res_G$sig <- if_else((res_G$log2FoldChange >1 & res_G$padj <0.01) | (res_G$log2FoldChange < -1 & res_G$padj <0.01), "TRUE", "FALSE")

res_H$sig <- if_else((res_H$log2FoldChange >1 & res_H$padj <0.01) | (res_H$log2FoldChange < -1 & res_H$padj <0.01), "TRUE", "FALSE")

res_I$sig <- if_else((res_I$log2FoldChange >1 & res_I$padj <0.01) | (res_I$log2FoldChange < -1 & res_I$padj <0.01), "TRUE", "FALSE")

res_J$sig <- if_else((res_J$log2FoldChange >1 & res_J$padj <0.01) | (res_J$log2FoldChange < -1 & res_J$padj <0.01), "TRUE", "FALSE")

res_K$sig <- if_else((res_K$log2FoldChange >1 & res_K$padj <0.01) | (res_K$log2FoldChange < -1 & res_K$padj <0.01), "TRUE", "FALSE")

# padj<0.01, Log2FC臒l?ݒ??ς?
#arrage by log2FC
graph_2 <- res_B %>% filter(sig == "TRUE")%>% dplyr::arrange(by = log2FoldChange)
graph_3 <- res_C %>% filter(sig == "TRUE")%>% dplyr::arrange(by = log2FoldChange)
graph_4 <- res_D %>% filter(sig == "TRUE")%>% dplyr::arrange(by = log2FoldChange)
graph_5 <- res_E %>% filter(sig == "TRUE")%>% dplyr::arrange(by = log2FoldChange)
graph_6 <- res_F %>% filter(sig == "TRUE")%>% dplyr::arrange(by = log2FoldChange)
graph_7 <- res_G %>% filter(sig == "TRUE")%>% dplyr::arrange(by = log2FoldChange)
graph_8 <- res_H %>% filter(sig == "TRUE")%>% dplyr::arrange(by = log2FoldChange)
graph_9 <- res_I %>% filter(sig == "TRUE")%>% dplyr::arrange(by = log2FoldChange)
graph_10 <- res_J %>% filter(sig == "TRUE")%>% dplyr::arrange(by = log2FoldChange)
graph_11 <- res_K %>% filter(sig == "TRUE")%>% dplyr::arrange(by = log2FoldChange)

#regulation 
graph_2_up <- dplyr::filter(graph_2, log2FoldChange>0) 
graph_2_down <- dplyr::filter(graph_2, log2FoldChange<0)
graph_3_up <- dplyr::filter(graph_3, log2FoldChange>0) 
graph_3_down <- dplyr::filter(graph_3, log2FoldChange<0) 
graph_4_up <- dplyr::filter(graph_4, log2FoldChange>0) 
graph_4_down <- dplyr::filter(graph_4, log2FoldChange<0) 
graph_5_up <- dplyr::filter(graph_5, log2FoldChange>0) 
graph_5_down <- dplyr::filter(graph_5, log2FoldChange<0)
graph_6_up <- dplyr::filter(graph_6, log2FoldChange>0) 
graph_6_down <- dplyr::filter(graph_6, log2FoldChange<0) 
graph_7_up <- dplyr::filter(graph_7, log2FoldChange>0) 
graph_7_down <- dplyr::filter(graph_7, log2FoldChange<0) 
graph_8_up <- dplyr::filter(graph_8, log2FoldChange>0) 
graph_8_down <- dplyr::filter(graph_8, log2FoldChange<0) 
graph_9_up <- dplyr::filter(graph_9, log2FoldChange>0) 
graph_9_down <- dplyr::filter(graph_9, log2FoldChange<0) 
graph_10_up <- dplyr::filter(graph_10, log2FoldChange>0) 
graph_10_down <- dplyr::filter(graph_10, log2FoldChange<0) 
graph_11_up <- dplyr::filter(graph_11, log2FoldChange>0) 
graph_11_down <- dplyr::filter(graph_11, log2FoldChange<0) 

# Count genes
# merge
downregulated <- dplyr::full_join(graph_11_down,graph_2_down) %>% full_join(graph_3_down) %>%
  full_join(graph_4_down)%>% full_join(graph_5_down)%>% full_join(graph_6_down)%>% full_join(graph_7_down)%>%
  full_join(graph_8_down)%>% full_join(graph_9_down)%>% full_join(graph_10_down)
dim(downregulated)
downregulated <- dplyr::distinct(downregulated,ENSEMBL)

upregulated <- dplyr::full_join(graph_11_up,graph_2_up) %>% full_join(graph_3_up) %>%
  full_join(graph_4_up)%>% full_join(graph_5_up)%>% full_join(graph_6_up)%>% full_join(graph_7_up)%>%
  full_join(graph_8_up)%>% full_join(graph_9_up)%>% full_join(graph_10_up)
dim(upregulated)
upregulated <- dplyr::distinct(upregulated,ENSEMBL)

#save
#Where: k:/haga/skin_vitro/data_analysis_R/DESeq2/
write.csv(downregulated, file = "20201030_HFF_invitro_downregulated_genes.csv")
write.csv(upregulated, file = "20201030_HFF_invitro_upregulated_genes.csv")

#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#intersect between in vivo and in vitro DEGs
# Set were to save data
setwd("D:/onedrive/大阪大学/論文/Figure/") 

# Check were you are now
getwd() 

# clear the decks
rm(list = ls())

# read(in vivo data)
cluster_1_2_down <- read.csv("D:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/DESeq2/20221209_cluster1_2_downregulated_genes.csv",header = TRUE, stringsAsFactors = FALSE,row.names = 1)
cluster_1_2_up <- read.csv("D:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/DESeq2/20221209_cluster1_2_upregulated_genes.csv",header = TRUE, stringsAsFactors = FALSE,row.names = 1)
cluster_1_3_down <- read.csv("D:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/DESeq2/20221209_cluster1_3_downregulated_genes.csv",header = TRUE, stringsAsFactors = FALSE,row.names = 1)
cluster_1_3_up <- read.csv("D:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/DESeq2/20221209_cluster1_3_upregulated_genes.csv",header = TRUE, stringsAsFactors = FALSE,row.names = 1)
cluster_2_3_down <- read.csv("D:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/DESeq2/20221209_cluster2_3_downregulated_genes.csv",header = TRUE, stringsAsFactors = FALSE,row.names = 1)
cluster_2_3_up <- read.csv("D:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/DESeq2/20221209_cluster2_3_upregulated_genes.csv",header = TRUE, stringsAsFactors = FALSE,row.names = 1)

# merge
downregulated_vivo <- dplyr::full_join(cluster_1_2_down,cluster_1_3_down) %>% full_join(cluster_2_3_down)
upregulated_vivo <- dplyr::full_join(cluster_1_2_up,cluster_1_3_up) %>% full_join(cluster_2_3_up)

# read(in vitro)
downregulated <- read.csv("k:/haga/skin_vitro/data_analysis_R/DESeq2/20201030_HFF_invitro_downregulated_genes.csv",header = TRUE, stringsAsFactors = FALSE,row.names = 1)
upregulated <- read.csv("k:/haga/skin_vitro/data_analysis_R/DESeq2/20201030_HFF_invitro_upregulated_genes.csv",header = TRUE, stringsAsFactors = FALSE,row.names = 1)

# vector
downregulated_vivo_vector <- downregulated_vivo$ENSEMBL
upregulated_vivo_vector <- upregulated_vivo$ENSEMBL
downregulated_vector <- downregulated$ENSEMBL
upregulated_vector <- upregulated$ENSEMBL

#Venn figure
#Downregulated
#Venn
color_palette <- c("#4F81BD","#C00000")
list <- list(in_vivo = downregulated_vivo_vector, in_vitro = downregulated_vector)
venn.diagram(list,
             filename="Fig1_H_Venn_intersect_Downregulated_vivo_vitro.png",
             category.names = c(expression(paste(italic("in vivo")), italic("in vitro"))),
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
#Upregulated

list <- list(in_vivo = upregulated_vivo_vector, in_vitro = upregulated_vector)
venn.diagram(list,
             filename="Fig1_H_Venn_intersect_Upregulated_vivo_vitro.png",
             category.names = c(expression(paste(italic("in vivo")), italic("in vitro"))),
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

# Common downreulated genes
common_down <- intersect(downregulated_vivo_vector ,downregulated_vector)

# Common upreulated genes
common_up <- intersect(upregulated_vivo_vector ,upregulated_vector)

#Count
# > sum(table(common_down))
# [1] 592
# > sum(table(common_up))
# [1] 502
sum(table(common_down))
sum(table(common_up))

# ENSEMBL to ENTREZID
down_ENTREZ <- bitr(common_down,  fromType = "ENSEMBL",
                    toType = c("SYMBOL", "ENTREZID"),
                    OrgDb = org.Hs.eg.db)

up_ENTREZ <- bitr(common_up,  fromType = "ENSEMBL",
                  toType = c("SYMBOL", "ENTREZID"),
                  OrgDb = org.Hs.eg.db)
#Save gene list
write.csv(down_ENTREZ, file = "Raw_data/Fig1/new_Downregulated_in_vivo_in_vitro_intersected_genes.csv")
write.csv(up_ENTREZ, file = "Raw_data/Fig1/new_Upregulated_in_vivo_in_vitro_intersected_genes.csv")
# clear the decks
rm(list = ls())

#Import intersected gene sets
down_ENTREZ <- read.csv("d:/onedrive/大阪大学/論文/Figure/Raw_data/Fig1/new_Downregulated_in_vivo_in_vitro_intersected_genes.csv",header = TRUE, stringsAsFactors = FALSE,row.names = 1)
up_ENTREZ <- read.csv("d:/onedrive/大阪大学/論文/Figure/Raw_data/Fig1/new_Upregulated_in_vivo_in_vitro_intersected_genes.csv",header = TRUE, stringsAsFactors = FALSE,row.names = 1)

#cluster
cluster <- list(Downregulated = down_ENTREZ$ENTREZID,Upregulated = up_ENTREZ$ENTREZID)

kegg <- compareCluster(geneCluster = cluster, fun = "enrichKEGG",  organism="hsa", 
                       keyType = "kegg",pAdjustMethod = "BH",minGSSize = 10, maxGSSize = 500,
                       pvalueCutoff  = 0.05,qvalueCutoff  = 0.05, use_internal_data = FALSE)
down_kegg <- enrichKEGG(gene = down_ENTREZ$ENTREZID,
                        organism = "hsa",
                        pAdjustMethod = "BH",
                        minGSSize = 10, maxGSSize = 500,
                        pvalueCutoff  = 0.05,qvalueCutoff  = 0.05, use_internal_data = FALSE)
up_kegg <- enrichKEGG(gene = up_ENTREZ$ENTREZID,
                      organism = "hsa",
                      pAdjustMethod = "BH",
                      minGSSize = 10, maxGSSize = 500,
                      pvalueCutoff  = 0.05,qvalueCutoff  = 0.05, use_internal_data = FALSE)
up_kegg
# Result 
dotplot(kegg, title ="KEGG analysis", font.size = 10)
write.csv(kegg, file ="Raw_data/Fig1/new_KEGG_HFF_compare.csv" )
write.csv(down_kegg, file ="Raw_data/Fig1/new_KEGG_HFF_down.csv" )
write.csv(up_kegg, file ="Raw_data/Fig1/new_KEGG_HFF_up.csv" )

#read
down <- read.csv("Raw_data/Fig1/KEGG_HFF_down.csv",header = T,stringsAsFactors = F) 
up <- read.csv("Raw_data/Fig1/KEGG_HFF_up.csv",header = T,stringsAsFactors = F)
down <- -log(down[,c(3,7)] %>% column_to_rownames(var="Description"))
up <- -log(up[,c(3,7)] %>% column_to_rownames(var="Description"))
cluster <- bind_rows(up,down)
cluster <- transform(cluster,regulation =c(rep("upregulated",nrow(up)),rep("downregulated",nrow(down))))
cluster

#Mark value over 12
cluster$sig <- if_else((cluster$p.adjust >= 12), "TRUE", "FALSE")
cluster

#save
write.csv(cluster,file = "Raw_data/Fig1/new_Data_for_KEGG_heatmap.csv")

# read
#Selected top 5 pathway for each regulation
cluster <- read.csv("Raw_data/Fig1/new_Data_for_KEGG_heatmap_changed.csv",header = T,stringsAsFactors = F,row.names = 1)

# pheatmap
BreaksList <- seq(0,12, by =3)
Col.Pal <- colorRampPalette(c("whitesmoke", "indianred"))
Col.List <- Col.Pal(length(BreaksList))

#Save
g1 <- pheatmap(mat = cluster,breaks = BreaksList,
               fontsize=7, fontsize_row = 6.5, fontsize_col = 6.5, fontfamily = "sans",
               cutree_cols =  1,cutree_rows = 1,
               color=Col.List,cluster_rows = F,cluster_cols = F,
               legend=F, legend_breaks=c(0,12),legend_labels =c("0","12 or greater"),
               main = "", angle_col = 45,
               treeheight_col = 0,  border_color = NA)
g1
graph2ppt(x=g1, file="Fig1_H_KEGG_heatmap", width = 2.3, height = 2)

#------------------------------------------------------------------------------------------------------------
# read
KEGG <- read.csv("Raw_data/Fig1/KEGG_HFF_compare.csv",header = TRUE, stringsAsFactors = FALSE, row.names = 1)
KEGG
KEGG <-filter(KEGG, ID =="hsa04350"|ID=="hsa04933")
KEGG_ID <-str_split(KEGG$geneID,"/")


as.data.frame(KEGG_ID[1])

A <- as.data.frame(KEGG_ID[1]) %>% dplyr::rename(ENTREZID =  c..1634....7057....7046....1030....7043....2331....92....64388...)
B <- as.data.frame(KEGG_ID[2]) %>% dplyr::rename(ENTREZID =  c..1958....1906....7046....5292....7043....3383....3569....6347...)

KEGG_ID <- full_join(A, B, by ="ENTREZID")
KEGG_ID <- distinct(KEGG_ID)

# ID convert
KEGG_genes <- bitr(KEGG_ID$ENTREZID,  fromType = "ENTREZID",
                   toType = c("SYMBOL" ,"ENSEMBL"),
                   OrgDb = org.Hs.eg.db)

SYMBOL <- KEGG_genes %>% dplyr::distinct(SYMBOL)
KEGG_genes
SYMBOL
# write 
write.csv(SYMBOL, file = "Raw_data/KEGG_fibrosis_genes_annotated.csv")

#------------------------------------------------------------------------------------------------------------
# Set were to save data
setwd("d:/onedrive/大阪大学/論文/Figure/") 

# clear the decks
rm(list = ls())

# read peak files
norm_count <- read.csv("D:/haga/human_skin/counts/clustering/TPM_genes_cuttoff_5/PCA/20200602_TPM_outlier_omitted.csv",header = TRUE,stringsAsFactors = FALSE) %>% dplyr::rename(ENSEMBL=X)
norm_count

# read in vivo, invitro fibrosis anotated genes
comon_fibrosis <- read.csv("Raw_data/KEGG_fibrosis_genes_annotated.csv",header = TRUE,stringsAsFactors = FALSE,row.names = 1) %>% as.data.frame()
fibrosis_genes <- bitr(comon_fibrosis$SYMBOL, fromType = "SYMBOL",
                       toType = c("ENSEMBL"),
                       OrgDb = org.Hs.eg.db)
fibrosis_genes

#join
TPM_fibrosis <- inner_join(fibrosis_genes,norm_count,by="ENSEMBL")%>%dplyr::select(-ENSEMBL)%>%column_to_rownames(var ="SYMBOL")
TPM_fibrosis

# Import age
age_info <- read.csv("D:/haga/human_skin/counts/clustering/old data/TPM_coding_genes_cutoff_5/cluster/10s_90s/20200415_Cluster_dendrogram_k=5_meancutoff_5_proteincoding.csv", header = TRUE, stringsAsFactors = FALSE)

age_info <- age_info %>% dplyr::select(-result)
age_info

# Add age data
TPM_fibrosis <-  t(TPM_fibrosis) %>% as.data.frame()
TPM_fibrosis <- rownames_to_column(TPM_fibrosis, var = "X")
TPM_fibrosis <- inner_join(age_info,TPM_fibrosis, by = "X")
TPM_fibrosis <- dplyr::arrange(TPM_fibrosis,age)
TPM_fibrosis
#set level
sample_factor <- c("young","young","young",
                   "young","young","young","young","young","young","young","young","young","young",
                   "young","young","young","young",
                   "middle","middle","middle","middle","middle","middle","middle","middle","middle","middle","middle","middle",
                   "middle","middle",
                   "middle","middle","middle","middle","middle","middle","middle","middle","middle","middle","middle","middle",
                   "old","old",
                   "old","old","old","old","old","old","old","old","old","old",
                   "old","old","old","old","old","old","old","old","old","old",
                   "old","old","old","old","old","old","old")
sample_factor2 <- factor(sample_factor, levels = c("young","middle","old"))

#Make group ID
TPM_fibrosis <- transform(TPM_fibrosis,group=sample_factor2)
TPM_fibrosis

# 
TPM <- TPM_fibrosis %>% column_to_rownames("X") %>% dplyr::select(-age, -group)
TPM
comon_fibrosis
comon_fibrosis$SYMBOL
mean <- data.frame(SYMBOL = comon_fibrosis$SYMBOL, mean =c(mean(TPM$DCN),mean(TPM$THBS1),mean(TPM$TGFBR1),mean(TPM$CDKN2B),mean(TPM$TGFB3),
                                                           mean(TPM$FMOD),mean(TPM$ACVR2A),mean(TPM$GREM2),mean(TPM$GDF6),
                                                           mean(TPM$EGR1),mean(TPM$EDN1),mean(TPM$PIM1),mean(TPM$ICAM1),
                                                           mean(TPM$IL6),mean(TPM$CCL2),mean(TPM$PLCE1)))
mean
# modify data
TPM_fibrosis <- TPM_fibrosis[,c(3:ncol(TPM_fibrosis)-1)]

# Calculate corr
rho <- cor(TPM_fibrosis, method = "spearman") %>% as.data.frame()
rho_age <- dplyr::select(rho,age)
rho_age <- dplyr::slice(rho_age, 2:nrow(rho_age))
rho_age$rank1 <- rank(-(rho_age$age))
rho_age <- rownames_to_column(rho_age, var ="SYMBOL")
rho_age <- arrange(rho_age,rank1)

rho_age
#merge
res <- dplyr::inner_join(rho_age,mean, by ="SYMBOL")
res$rank2 <- rank(-(res$mean))
res <- mutate(res, rank = (rank1+rank2)/2)
res$rank_res <- rank(res$rank)                  
res <- arrange(res,rank_res)
res <- res[,c(1,2,4,7)]
res <- dplyr::rename(res, spearman =age)
res <- mutate(res, normalize_LogTPM = log2(mean+1))
res
#save
write.csv(res, file = "Raw_data/in_vivo_genes_fibrosis.csv")

# Mark Top rank genes
res$Top5_genes <- if_else(res$rank_res<6, "TRUE", "FALSE")
res

#in vitro 
norm_count_vitro <- read.csv("k:/haga/skin_vitro/data_analysis_R/TPM/20201030_TPM_normalized_data_cutoff_mean_5.csv",header = TRUE,stringsAsFactors = FALSE) %>% dplyr::rename(ENSEMBL=X)
head(norm_count_vitro)

#inner join
merged <- inner_join(fibrosis_genes,norm_count_vitro, by="ENSEMBL")
HFF_sample <- merged[,c(1,11:ncol(merged))]
HFF_sample
t_data <- HFF_sample %>% column_to_rownames(var="SYMBOL") %>% t() %>%as.data.frame() %>% rownames_to_column(var="sample")
t_data
#PDL info
PDL <- data.frame(sample=t_data$sample,PDL =c(16,16,16,26,26,26,46,46,46,64,64,64,74,74,74))
res_temp <- inner_join(PDL,t_data)
res_temp <- res_temp[,c(2:ncol(res_temp))]
rho_vitro <- cor(res_temp, method = "spearman") %>% as.data.frame()
rho_PDL <- dplyr::select(rho_vitro,PDL)
rho_PDL <- dplyr::slice(rho_PDL, 2:nrow(rho_PDL))
rho_PDL <- rho_PDL %>% rownames_to_column(var="SYMBOL")
rho_PDL

#merge with vivo result
res <- inner_join(rho_PDL,res, by="SYMBOL")
res

#Rank by correlation
res <- res[,c(1:3)]
res$rank1 <-rank(-(res$PDL))
res$rank2 <-rank(-(res$spearman))
res <- mutate(res, rank = (rank1+rank2)/2)
res$res <- rank((res$rank))
res <- arrange(res,by=res)
res
# Mark Top rank genes
res$genes <- if_else(res$SYMBOL=="THBS1"|res$SYMBOL=="FMOD", "TRUE", "FALSE")

#make data
write.csv(res,file = "Raw_data/Fig1/Correlation_data_in_vivo_in_vitro.csv")

# clear the decks
rm(list = ls())

#Import data
res <- read.csv("Raw_data/Fig1/Correlation_data_in_vivo_in_vitro.csv",header = T,stringsAsFactors = F,row.names = 1)
add_data <- data.frame("Group"=c(rep("TGF",nrow(res))))
res <- cbind(res,add_data)

#Plot
g<- ggplot(res, aes(x = spearman, y = PDL, label =SYMBOL, color = SYMBOL,group=Group)) 
g<- g + geom_point(size =0.5, alpha =0.5) 
g<- g + geom_smooth(method="lm",formula='y~x', colour = "gray75", se=FALSE, size = 0.5,alpha=0.5)
g<- g + stat_poly_eq(formula = y ~ x,
                     aes(label = paste(stat(eq.label),
                                       stat(rr.label),
                                       sep = "~~~")),
                     parse = TRUE)
g<- g + theme_bw() 
g<- g + xlab (expression(paste(italic("in vivo "), "(age)"))) + ylab(expression(paste(italic("in vitro "), "(PDL)"))) 
g <- g + geom_text_repel(fontent.size = 1.5,show.legend = FALSE,size=1.5) 
g <- g + scale_x_continuous(expand = c(0, 0), limits = c(-0.2,0.7)) 
g <- g + scale_y_continuous(expand = c(0, 0), limits = c(-0.2,1.2)) 
g <- g + scale_colour_manual(values=c("THBS1"="#C00000",
                                      "FMOD"="#4F81BD",
                                      "TGFBR1"="gray50",
                                      "CDKN2B"="gray50",
                                      "DCN"="gray50",
                                      "GDF6"="gray50",
                                      "PLCE1"="gray50",
                                      "GREM2"="gray50",
                                      "TGFB3"="gray50",
                                      "ICAM1"="gray50",
                                      "IL6"="gray50",
                                      "CCL2"="gray50",
                                      "PIM1"="gray50",
                                      "EDN1"="gray50",
                                      "EGR1"="gray50",
                                      "ACVR2A"="gray50")) 

g <- g + guides(colour=FALSE) 
g <- g + theme(axis.title = element_text(size = 8), 
               axis.text = element_text(size = 8), 
               axis.text.x = element_text(size = 8), 
               axis.text.y = element_text(size = 8), 
               legend.text = element_text(size = 8), 
               legend.title = element_text(size = 8)) + theme(plot.title = element_text(size = 8))

g
graph2ppt(x=g, file="Fig1_J_correlation_vivo_vitro", width = 2.1, height = 2.1)
