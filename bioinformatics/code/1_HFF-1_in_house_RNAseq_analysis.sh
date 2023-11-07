#!/usr/bin/bash

echo "####################################################################"
echo "#              Dataset for Dermal fibroblast senescence            #"
echo "#              Analysis of HFF1                                    #"
echo "#              Data prepared by Haga et al.                        #"
echo "#              Data type:    RNA-seq (Paired-end)                  #"
echo "#              PDL: PDL24, 36, 47                                  #"
echo "#              Treatment: Control, TGF-beta1                       #"
echo "#              Data availability: PRJDB15707                       #"
echo "#              N=3                                                 #"
echo "####################################################################"

PROJECT_PATH=/mnt/d/haga/HFF_RS_TGF
#Current directory (Directory containing fastq files)
#/mnt/d/haga/HFF_RS_TGF/fastq

for SAMPLE in `ls *_1.fq.gz`
do
    name=${SAMPLE%_1.fq.gz}
#trim_galore fastqc
    echo "Perform trim galore(version 0.6.6)"
    trim_galore --fastqc --cores 30 --paired ${PROJECT_PATH}/fastq/${name}_1.fq.gz ${PROJECT_PATH}/fastq/${name}_2.fq.gz \
    -o ${PROJECT_PATH}/trimmed_fastq

#Mapping using hisat2
    echo "Perform hisat2(version 2.2.1)"
    hisat2 --summary-file ${PROJECT_PATH}/hisat2/${name}.txt -p 30 -q \
    -x ${PROJECT_PATH}/hisat2/grch38/genome \
    -1 ${PROJECT_PATH}/trimmed_fastq/${name}_1_val_1.fq.gz \
    -2 ${PROJECT_PATH}/trimmed_fastq/${name}_2_val_2.fq.gz \
    -S ${PROJECT_PATH}/SAM/${name}.sam
    
#samtools
    echo "Perform samtools(ver. 1.9)" 
    samtools view -@ 30 -b -q 10 ${PROJECT_PATH}/SAM/${name}.sam > ${PROJECT_PATH}/BAM/${name}.bam
    samtools sort -@ 30 ${PROJECT_PATH}/BAM/${name}.bam > ${PROJECT_PATH}/BAM/${name}_sorted.bam
    samtools index -@ 30 ${PROJECT_PATH}/BAM/${name}_sorted.bam

#featureCounts
cd /mnt/d/haga/HFF_RS_TGF/BAM
echo "Perform featureCounts(subread ver. 2.0.1)"
output=20211025_HFF_counts.txt
    featureCounts -p -T 30 -t exon -g gene_id \
    -a ${PROJECT_PATH}/index/gencode.v33.annotation.gtf \
    -o ${PROJECT_PATH}/featureCounts/${output} *_sorted.bam

done