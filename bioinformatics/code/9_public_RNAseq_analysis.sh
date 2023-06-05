#!/usr/bin/bash

echo "####################################################################"
echo "#              Dataset for in vitro dermal fibroblast senescence   #"
echo "#              Analysis of pubilic HFF1                            #"
echo "#              citation: Marthandan et al., 2016                   #"
echo "#              Data type:    RNA-seq (single-end)                  #"
echo "#              PDL: PDL16, 26, 46, 64, 74                          #"
echo "#              N=3                                                 #"
echo "####################################################################"

<< COMMENTOUT
Get SRA files from DBCL
study: GSE63577
secondary accession: SRP050179
COMMENTOUT
#HFF_PD16
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX766/SRX766250/SRR1660543/SRR1660543.sra
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX766/SRX766251/SRR1660544/SRR1660544.sra
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX766/SRX766252/SRR1660545/SRR1660545.sra
#HFF_PD26
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX136/SRX1362168/SRR2751110/SRR2751110.sra
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX136/SRX1362169/SRR2751111/SRR2751111.sra
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX136/SRX1362170/SRR2751112/SRR2751112.sra
#HFF_PD46
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX136/SRX1362171/SRR2751113/SRR2751113.sra
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX136/SRX1362172/SRR2751114/SRR2751114.sra
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX136/SRX1362173/SRR2751115/SRR2751115.sra
# HFF_PD64
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX136/SRX1362174/SRR2751116/SRR2751116.sra
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX136/SRX1362176/SRR2751117/SRR2751117.sra
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX136/SRX1362177/SRR2751118/SRR2751118.sra
# HFF_PD74
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX766/SRX766253/SRR1660546/SRR1660546.sra
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX766/SRX766254/SRR1660547/SRR1660547.sra
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX766/SRX766255/SRR1660548/SRR1660548.sra

SAMPLES="SRR1660543 SRR1660544 SRR1660545 \
SRR1660546 SRR1660548 SRR1660547 SRR2751110 SRR2751111 SRR2751112 SRR2751113 SRR2751114 SRR2751115 SRR2751116 SRR2751117 SRR2751118"
for SAMPLE in $SAMPLES; 
do
#fastq-dump (sra to fastq)
fastq-dump /mnt/d/haga/invitro/sample/${SAMPLE}.sra
#trim_galore fastqc
trim_galore --fastqc --cores 4 ./${SAMPLE}.fastq -o ../trimmed_fastq
#Mapping using hisat2
hisat2 --summary-file ./${SAMPLE}.txt -p 16 -q -x ./genome_index/genome_index -U ../trimmed_fastq/${SAMPLE}_trimmed.fq -S /mnt/d/haga/invitro/hisat2/${SAMPLE}.sam
#samtools
samtools view -@ 16 -b -q 10 /mnt/d/haga/invitro/hisat2/${SAMPLE}.sam > /mnt/d/haga/invitro/BAM/${SAMPLE}.bam
samtools sort -@ 16 /mnt/d/haga/invitro/BAM/${SAMPLE}.bam > /mnt/d/haga/invitro/BAM/${SAMPLE}_sorted.bam
samtools index -@ 16 /mnt/d/haga/invitro/BAM/${SAMPLE}_sorted.bam
#featureCounts
featureCounts --verbose -T 16 -t exon -g gene_id -a ../index/gencode.v33.annotation.gtf -o ./output_counts.txt ../BAM/*_sorted.bam

done

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

echo "####################################################################"
echo "#              Dataset for in vivo dermal fibroblast senescence    #"
echo "#              Analysis of human skin fibroblast                   #"
echo "#              citation: Fleischer et al., 2018                    #"
echo "#              Data type:    RNA-seq (single-end)                  #"
echo "####################################################################"

SAMPLES="SRR7093809 SRR7093810 SRR7093811 SRR7093812 SRR7093813 \
SRR7093814 SRR7093815 SRR7093816 SRR7093817 SRR7093818 \
SRR7093819 SRR7093820 SRR7093821 SRR7093822 SRR7093823 \
SRR7093824 SRR7093825 SRR7093826 SRR7093827 SRR7093828 \
SRR7093829 SRR7093830 SRR7093831 SRR7093832 SRR7093833 \
SRR7093834 SRR7093835 SRR7093836 SRR7093837 SRR7093838 \
SRR7093839 SRR7093840 SRR7093841 SRR7093842 SRR7093843 \
SRR7093844 SRR7093845 SRR7093846 SRR7093847 SRR7093848 \
SRR7093849 SRR7093850 SRR7093851 SRR7093852 SRR7093853 \
SRR7093854 SRR7093855 SRR7093856 SRR7093857 SRR7093858 \
SRR7093859 SRR7093860 SRR7093861 SRR7093862 SRR7093863 \
SRR7093864 SRR7093865 SRR7093866 SRR7093867 SRR7093868 \
SRR7093869 SRR7093870 SRR7093871 SRR7093872 SRR7093873 \
SRR7093874 SRR7093875 SRR7093876 SRR7093877 SRR7093878 \
SRR7093879 SRR7093880 SRR7093881 SRR7093882 SRR7093883 \
SRR7093884 SRR7093885 SRR7093886 SRR7093887 SRR7093888 \
SRR7093889 SRR7093890 SRR7093891 SRR7093892 SRR7093893 \
SRR7093894 SRR7093895 SRR7093896 SRR7093897 SRR7093898 \
SRR7093899 SRR7093900 SRR7093901 SRR7093902 SRR7093903 \
SRR7093904 SRR7093905 SRR7093906 SRR7093907 SRR7093908 \
SRR7093909 SRR7093910 SRR7093911 SRR7093912 SRR7093913 \
SRR7093914 SRR7093915 SRR7093916 SRR7093917 SRR7093918 \
SRR7093919 SRR7093920 SRR7093921 SRR7093922 SRR7093923 \
SRR7093924 SRR7093925 SRR7093926 SRR7093927 SRR7093928 \
SRR7093929 SRR7093930 SRR7093931 SRR7093932 SRR7093933 \
SRR7093934 SRR7093935 SRR7093936 SRR7093937 SRR7093938 \
SRR7093939 SRR7093940 SRR7093941 SRR7093942 SRR7093943 \
SRR7093944 SRR7093945 SRR7093946 SRR7093947 SRR7093948 \
SRR7093949 SRR7093950 SRR7093951"
for SAMPLE in $SAMPLES; do
#Download sra files
prefetch ${SAMPLE}
#fastq-dump (sra to fastq)
fastq-dump ${SAMPLE}.sra
#trim_galore fastqc
trim_galore --fastqc ./SRA/${SAMPLE}.sra.fastq -o ../trimmed_fastq
#Mapping using hisat2
hisat2 -p 6 -q -x ../index/grch38/genome -U ../trimmed_fastq/${SAMPLE}.sra_trimmed.fq -S ./mapped_sam/${SAMPLE}.sam
#samtools
samtools view -S -b ./mapped_sam/${SAMPLE}_trimmed.sam > ./mapped_bam/${SAMPLE}.bam
samtools sort ./mapped_bam/${SAMPLE}.bam > ./mapped_bam/${SAMPLE}_sorted.bam
samtools index ./mapped_bam/${SAMPLE}_sorted.bam
#featureCounts
featureCounts --verbose -T 5 -t exon -g gene_id -a ./gene_annotation/gencode.v33.annotation.gtf -o ./counts/counts.txt ./mapped_bam/*_sorted.bam

done