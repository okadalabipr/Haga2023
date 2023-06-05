#!/usr/bin/bash

echo "####################################################################"
echo "#              Dataset for Dermal fibroblast senescence            #"
echo "#              Analysis of HFF1                                    #"
echo "#              Data prepared by Haga et al.                        #"
echo "#              Data type:    ATAC-seq, ChIP-seq (Paired-end)       #"
echo "#              PDL: PDL24, 36, 47                                  #"
echo "#              N=2                                                 #"
echo "####################################################################"

#Concatenated using “cat” command ATAC-seq differential peaks
#Use nf-core pipeline output
cd /mnt/d/HFF_RS_TGF_ATAC/results/bwa/mergedLibrary/macs/broadPeak/consensus/Diff_ATAC_PDL/
cat HFF_PDL_24_CONvsHFF_PDL_36_CON.mLb.clN.deseq2.FDR0.05.results.bed HFF_PDL_24_CONvsHFF_PDL_47_CON.mLb.clN.deseq2.FDR0.05.results.bed HFF_PDL_36_CONvsHFF_PDL_47_CON.mLb.clN.deseq2.FDR0.05.results.bed > tmp1.bed
sortBed -i tmp1.bed > tmp2.bed
mergeBed -i tmp2.bed > Differential_ATAC_output.bed
rm -rf tmp1.bed tmp2.bed

#Concatenated using “cat” command ChIP-seq(H3K27Ac) differential peaks
#Use nf-core pipeline output
cd HFF_RS_TGF_ChIP/results/bwa/mergedLibrary/macs/broadPeak/consensus/Diff_H3K27Ac_PDL/
cat PDL24_CON_AcvsPDL36_CON_Ac.deseq2.FDR0.05.results.bed PDL24_CON_AcvsPDL47_CON_Ac.deseq2.FDR0.05.results.bed PDL36_CON_AcvsPDL47_CON_Ac.deseq2.FDR0.05.results.bed > tmp1.bed
sortBed -i tmp1.bed > tmp2.bed
mergeBed -i tmp2.bed > Differential_H3K27Ac_output.bed
rm -rf tmp1.bed tmp2.bed
