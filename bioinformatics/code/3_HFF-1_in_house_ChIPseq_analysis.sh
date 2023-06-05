#!/usr/bin/bash

echo "####################################################################"
echo "#              Dataset for Dermal fibroblast senescence            #"
echo "#              Analysis of HFF1                                    #"
echo "#              Data prepared by Haga et al                         #"
echo "#              Data type:    ChIP-seq (Paired-end)                 #"
echo "#              PDL: PDL24, 36, 47                                  #"
echo "#              N=2                                                 #"
echo "####################################################################"

# for more info on settings: https://nf-co.re/chipseq/1.2.2
nextflow run nf-core/chipseq \
    -r 1.2.2 \
    -profile singularity \
    --input Samplesheet_ChIP.csv \
    --genome GRCh38 \
    --save_reference \
    --max_cpus 16 \
    --max_memory 256.GB
