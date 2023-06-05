# Bioinformatic analysis of skin aging
## Description
This repository contains the source code for the sequence analysis.
Please refer to files in [`ref_file/`](./ref_file/). Required files are listed in each script.

Primary analysis for RNA-seq was performed using script bellow.
* RNA-seq -> `1_HFF-1_in_house_RNAseq_analysis.sh`

Primary analysis for ATAC-seq and ChIP-seq were performed using an established pipeline.
* ATAC-seq -> https://nf-co.re/atacseq/1.2.1
* ChIP-seq -> https://nf-co.re/chipseq/1.2.2

## Requirements
### Software and Algorithms
Ubuntu version 20.04
| Package Name                                                   | Version                                                                                                  |
| -------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------- |
| Trim Galore!                                                   | 0.6.6                                                                                                    |
| hisat2        	                                             | 2.2.1                                                                                                    |
| Samtools                                                       | 1.9                                                                                                      |
| Subread         	                                             | 2.0.1                                                                                                    |
| BEDTools      	                                             | 2.30.0                                                                                                   |
| HOMER         	                                             | 4.11                                                                                                     |

R version 4.2.1
| Package Name                                                   | Version                                                                                                  |
| -------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------- |
| DoRothEA                                                       | 1.8.0                                                                                                    |
| clusterProfiler	                                             | 4.4.4                                                                                                    |
| DESeq2	                                                     | 1.36.0                                                                                                   |
| ChIPseeker	                                                 | 1.32.1                                                                                                   |
| rrcov	                                                         | 1.5.2                                                                                                    |

## Usage
To reproduce results, please perform the analysis in the following order.
* `1_HFF-1_in_house_RNAseq_analysis.sh`
* `2_HFF-1_in_house_ATACseq_analysis.sh`
* `3_HFF-1_in_house_ChIPseq_analysis.sh`
* `4_Transcription_factor_enrichment.R`
* `5_Gain_peaks_ATAC_ChIP_seq.R`
* `6_Motif_enrichment_ATAC.sh`
* `7_Differential_peak_ATAC_ChIP.sh`
* `8_Gene_annotation_differential_peak_ATAC_ChIP.R`
* `9_public_RNAseq_analysis.sh`
* `10_heatmap_of_KEGG_analysis_invivo_invitro.R`
