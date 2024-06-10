
##Convert to bam

#samtools sort -@ 5 bwa/$1/$1.sam -o bwa/$1/$1.bam -T bwa/$1/$1 -O BAM

## Get all mapped and unmapped reads, remove non-primary alignments
samtools view -F 256 -hb -o bwa/$1/$1.total.bam bwa/$1/$1.bam
samtools index bwa/$1/$1.total.bam
samtools flagstat bwa/$1/$1.total.bam > bwa/$1/$1.total.flagstat
fastqc bwa/$1/$1.total.bam
#rm bwa/$1/$1.sam
rm bwa/$1/$1.bam

##Get primary alignments -!!!PAIRED_END!!!
samtools view -Sbh -F 260 -f 3 bwa/$1/$1.total.bam -o bwa/$1/$1.primary.bam
samtools index bwa/$1/$1.primary.bam
samtools flagstat bwa/$1/$1.primary.bam > bwa/$1/$1.primary.flagstat
fastqc bwa/$1/$1.primary.bam

##filter for uniquely mapped reads (mapping quality > 20)
samtools view -b -q 20 bwa/$1/$1.primary.bam -o bwa/$1/$1.unique.bam
samtools index bwa/$1/$1.unique.bam
samtools flagstat bwa/$1/$1.unique.bam > bwa/$1/$1.unique.flagstat
fastqc bwa/$1/$1.unique.bam

##remove duplicates with picard from primary
picard-2.24.0 MarkDuplicates AS=true I=bwa/$1/$1.primary.bam O=bwa/$1/$1.primary_NR.bam M=bwa/$1/$1.primary_NR.dupmetrics REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT > bwa/$1/$1.dups.log
samtools index bwa/$1/$1.primary_NR.bam 2> bwa/$1/$1.dups.log
fastqc bwa/$1/$1.primary_NR.bam 2> bwa/$1/$1.dups.log
samtools flagstat bwa/$1/$1.primary_NR.bam > bwa/$1/$1.primary_NR.flagstat 2> bwa/$1/$1.dups.log

##remove duplicates with picard from uniq
picard-2.24.0 MarkDuplicates AS=true I=bwa/$1/$1.unique.bam O=bwa/$1/$1.unique_NR.bam M=bwa/$1/$1.unique_NR.dupmetrics REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT > bwa/$1/$1.dups.unique.log
samtools index bwa/$1/$1.unique_NR.bam 2> bwa/$1/$1.dups.unique.log
fastqc bwa/$1/$1.unique_NR.bam 2> bwa/$1/$1.dups.unique.log
samtools flagstat bwa/$1/$1.unique_NR.bam > bwa/$1/$1.unique_NR.flagstat 2> bwa/$1/$1.dups.unique.log

multiqc -f bwa/$1 -o bwa/$1






