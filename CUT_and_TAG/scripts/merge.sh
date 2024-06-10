
## Merge bam

samtools merge -@ 5 $1 -o $2.bam
samtools index $2.bam

