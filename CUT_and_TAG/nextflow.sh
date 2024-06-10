export NXF_OPTS='-Xms10g -Xmx40g'

## NF-core cut and run pipeline

nextflow run nf-core/cutandrun -profile singularity -r 3.1 -with-trace -with-tower --input samples.csv --genome GRCm38 --outdir results --save_reference --save_merged_fastq --save_trimmed --save_align_intermed --trim_nextseq 20 --minimum_alignment_q_score 0 --normalisation_mode CPM --normalisation_binsize 1 -c nextflow.config


