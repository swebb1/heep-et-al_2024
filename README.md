# heep_et_al_2024

## Two-factor authentication underpins the precision of piRNA-directed LINE1 DNA methylation

Source code for analysing ChIP-seq, ATAC-seq and CUT&TAG data.

ChIP-seq and ATAC-seq were processed with a custom Snakemake workflow.
CUT&TAG data was processed with NF-core CUT&RUN Nextflow pipeline v3.1

Additional scripts were used to remove duplicate reads, merge replicates, create log2 bigWig files and call peaks

Heatmaps over repetitive elements were generated with Deeptools.
Analysis of peaks was performed in R with custom scripts.
