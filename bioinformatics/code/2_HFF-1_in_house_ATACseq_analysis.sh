#!/usr/bin/bash

echo "####################################################################"
echo "#              Dataset for Dermal fibroblast senescence            #"
echo "#              Analysis of HFF1                                    #"
echo "#              Data prepared by Haga et al                         #"
echo "#              Data type:    ATAC-seq (Paired-end)                 #"
echo "#              PDL: PDL24, 36, 47                                  #"
echo "#              Treatment: Control, TGF-beta1(Only for PDL24)       #"
echo "#              Data availability: PRJDB15707                       #"
echo "#              N=2                                                 #"
echo "####################################################################"

# for more info on settings: https://nf-co.re/atacseq/1.2.1
nextflow run nf-core/atacseq \
    -r 1.2.1 \
    -profile singularity \
    --input Samplesheet_atac.csv \
    --genome GRCh38 \
    --save_reference \
    --max_cpus 8 \
    --max_memory 256.GB
