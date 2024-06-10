
## Call Peaks

macs2 callpeak -t $1 -c $2 -f BAMPE -g mm -n $3

macs2 callpeak -t $1 -c $2 -f BAMPE -g mm -n $3_all --keep-dup all

## Get read counts at peaks

multiBamSummary BED-file --BED $3_peaks.narrowPeak --bamfiles $1 -o $3_peaks_sum.npz --outRawCounts $3_peaks_sum.tsv -e -p 10

multiBamSummary BED-file --BED $3_all_peaks.narrowPeak --bamfiles $1 -o $3_peaks_sum.all.npz --outRawCounts $3_peaks_sum.all.tsv -e -p 10

