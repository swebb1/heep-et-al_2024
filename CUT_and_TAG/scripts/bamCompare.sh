
bamCoverage -b $1.bam -o $1.bw  -bs 1 -p 5 --normalizeUsing CPM --ignoreForNormalization MT
bamCoverage -b $1_NR.bam -o $1_NR.bw -bs 1 -p 5 --normalizeUsing CPM --ignoreForNormalization MT

bamCompare -b1 $1.bam -b2 $2.bam -o $1.log2.bw -bs 1 -p 5 --normalizeUsing BPM --exactScaling --ignoreForNormalization MT -e --scaleFactorsMethod None
bamCompare -b1 $1_NR.bam -b2 $2_NR.bam -o $1_NR.log2.bw -bs 1 -p 5 --normalizeUsing BPM --exactScaling --ignoreForNormalization MT -e --scaleFactorsMethod None

