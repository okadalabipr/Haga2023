# 5_Gain_peaks_ATAC_seq
# Obtain gained peaks with PDL increase for ATAC-seq and ChIP-seq (H3K27Ac)

# Libraies I need
library(tidyverse)

#Set were to save data
setwd("D:/haga/github/Haga2023/bioinformatics/ref_file/") 

# Check were you are now
getwd() 

# clear the decks
rm(list = ls())

# Read signal peaks
# read Differential peak files(use nf-core output)
peak_signal<-read.csv("D:/haga/HFF_RS_TGF_ATAC/results/bwa/mergedLibrary/macs/broadPeak/consensus/deseq2/consensus_peaks.mLb.clN.results.txt",header=T,sep="\t", stringsAsFactors = F)
colnames(peak_signal)

#PDL24_vs_PDL36
#adjp<0.05
#Gained(Log2FC>0) 22
PDL24_vs_PDL36 <- peak_signal[,c(1:6,13:18)]
PDL24_vs_PDL36 <- PDL24_vs_PDL36%>%filter(HFF_PDL_24_CONvsHFF_PDL_36_CON.padj<0.05)
PDL24_vs_PDL36_gain <- PDL24_vs_PDL36%>%filter(HFF_PDL_24_CONvsHFF_PDL_36_CON.log2FoldChange>0)

#Save
write.table(PDL24_vs_PDL36_gain[,2:4],"ATAC_PDL24_PDL36_gain.bed",quote=F,sep="\t",row.names=F,col.names=F)

#PDL24_vs_PDL47
#adjp<0.05
#Gained(Log2FC>0) 1029
PDL24_vs_PDL47 <- peak_signal[,c(1:6,25:30)]
PDL24_vs_PDL47 <- PDL24_vs_PDL47%>%filter(HFF_PDL_24_CONvsHFF_PDL_47_CON.padj<0.05)
PDL24_vs_PDL47_gain <- PDL24_vs_PDL47%>%filter(HFF_PDL_24_CONvsHFF_PDL_47_CON.log2FoldChange>0)

#Save
write.table(PDL24_vs_PDL47_gain[,2:4],"ATAC_PDL24_PDL47_gain.bed",quote=F,sep="\t",row.names=F,col.names=F)

#PDL36_vs_PDL47
#adjp<0.05
#Gained(Log2FC>0) 21
PDL36_vs_PDL47 <- peak_signal[,c(1:6,67:72)]
PDL36_vs_PDL47 <- PDL36_vs_PDL47%>%filter(HFF_PDL_36_CONvsHFF_PDL_47_CON.padj<0.05)
PDL36_vs_PDL47_gain <- PDL36_vs_PDL47%>%filter(HFF_PDL_36_CONvsHFF_PDL_47_CON.log2FoldChange>0)
dim(PDL36_vs_PDL47_gain)

#Save
write.table(PDL36_vs_PDL47_gain[,2:4],"ATAC_PDL36_PDL47_gain.bed",quote=F,sep="\t",row.names=F,col.names=F)

#------------------------------------------------------------------------------------------------------------
# Obtain DESEQ2 output and make bed files for motif enrichment
# ChIP-seq (H3K27Ac)

#Set were to save data
setwd("D:/haga/github/Haga2023/bioinformatics/ref_file/") 

# Check were you are now
getwd() 

# clear the decks
rm(list = ls())

#Read signal peaks
peak_signal<-read.csv("D:/haga/HFF_RS_TGF_ChIP/results/bwa/mergedLibrary/macs/broadPeak/consensus/H3K27Ac/deseq2/H3K27Ac.consensus_peaks.results.txt",header=T,sep="\t", stringsAsFactors = F)
colnames(peak_signal)

#PDL24_vs_PDL36
#adjp<0.05
#Gained(Log2FC>0) 1030
PDL24_vs_PDL36 <- peak_signal[,c(1:6,13:18)]
PDL24_vs_PDL36 <- PDL24_vs_PDL36%>%filter(PDL24_CON_AcvsPDL36_CON_Ac.padj<0.05)
PDL24_vs_PDL36_gain <- PDL24_vs_PDL36%>%filter(PDL24_CON_AcvsPDL36_CON_Ac.log2FoldChange>0)

#Save
write.table(PDL24_vs_PDL36_gain[,2:4],"H3K27Ac_PDL24_PDL36_gain.bed",quote=F,sep="\t",row.names=F,col.names=F)

#PDL24_vs_PDL47
#adjp<0.05
#Gained(Log2FC>0) 925
PDL24_vs_PDL47 <- peak_signal[,c(1:6,25:30)]
PDL24_vs_PDL47 <- PDL24_vs_PDL47%>%filter(PDL24_CON_AcvsPDL47_CON_Ac.padj<0.05)
PDL24_vs_PDL47_gain <- PDL24_vs_PDL47%>%filter(PDL24_CON_AcvsPDL47_CON_Ac.log2FoldChange>0)

#Save
write.table(PDL24_vs_PDL47_gain[,2:4],"H3K27Ac_PDL24_PDL47_gain.bed",quote=F,sep="\t",row.names=F,col.names=F)

#PDL36_vs_PDL47
#adjp<0.05
#Gained(Log2FC>0) 31
PDL36_vs_PDL47 <- peak_signal[,c(1:6,67:72)]
PDL36_vs_PDL47 <- PDL36_vs_PDL47%>%filter(PDL36_CON_AcvsPDL47_CON_Ac.padj<0.05)
PDL36_vs_PDL47_gain <- PDL36_vs_PDL47%>%filter(PDL36_CON_AcvsPDL47_CON_Ac.log2FoldChange>0)

#Save
write.table(PDL36_vs_PDL47_gain[,2:4],"H3K27Ac_PDL36_PDL47_gain.bed",quote=F,sep="\t",row.names=F,col.names=F)
