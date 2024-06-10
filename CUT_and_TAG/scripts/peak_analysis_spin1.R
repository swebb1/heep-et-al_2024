library(GenomicRanges)
library(genomation)
library(ChIPseeker)
library(stringr)
library(dplyr)
library(readr)
library(purrr)
library(ggplot2)

## Read in combined peaks
s1f <- list.files("../../bam_files/CUT_and_TAG_peaks",pattern = "SPIN1_E145_NR_peaks.narrowPeak",full.names = T,recursive = F)
names(s1f) <- s1f |> str_remove("../../bam_files/CUT_and_TAG_peaks/") |> str_remove("_peaks.narrowPeak")

s1.gr <- readPeakFile(s1f)

## Read in individual peaks
s1f <- list.files("../../bam_files/CUT_and_TAG_peaks",pattern = "SPIN1_E145_R[1-3]_NR_peaks.narrowPeak",full.names = T,recursive = F)
names(s1f) <- s1f |> str_remove("../../bam_files/CUT_and_TAG_peaks/") |> str_remove("_peaks.narrowPeak")

s1.grl <- s1f |> map(~readPeakFile(.x))
names(s1.grl) <- names(s1f)

# Subset peaks by score, enrichment and no of reads
pt = 191.8 # Peak threshold
rct = 20

s1rc <- read_tsv("../../bam_files/CUT_and_TAG_peaks/SPIN1_E145_peaks_sum.NR.tsv",col_names = c("chr","start","end","SPIN1_E145_R1","SPIN1_E145_R2","SPIN1_E145_R3"),skip = 1) |>
  mutate(start=start+1,length=end-start)

rcs1 = left_join(as.data.frame(s1.gr),s1rc,by = join_by(seqnames==chr,start==start,end==end)) |>
  filter(SPIN1_E145_R1 >= rct) |> pull(V4)

s1.f.gr <- subset(s1.gr,s1.gr$V5 >= pt & s1.gr$V4 %in% rcs1 & width(s1.gr) < 8000)

s1.f.gr$hits = 0
for(i in 1:length(s1.grl)){
  qh = findOverlaps(s1.f.gr,s1.grl[[i]]) |> queryHits()
  s1.f.gr$hits[qh] = s1.f.gr$hits[qh] + 1 
}
s1.f.gr <- subset(s1.f.gr,s1.f.gr$hits>=2)

seqlevelsStyle(s1.gr) <- "UCSC"
seqlevelsStyle(s1.f.gr) <- "UCSC"
saveRDS(s1.gr,"s1.peaks.rds")
saveRDS(s1.f.gr,"s1.peaks.f.rds")

## Annotate regions
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

txdb<-TxDb.Mmusculus.UCSC.mm10.knownGene
ens.txdb<-makeTxDbFromGFF("/datastore/home/genomes/mouse/mm10/annotation/Mus_musculus.GRCm38.78.edited.gtf",format = "gtf",organism = "Mus musculus")

peakAnno_s1<-annotatePeak(peak = s1.gr, TxDb=ens.txdb,
                          tssRegion=c(-3000, 3000),level="gene")

saveRDS(peakAnno_s1,"peakAnno_s1.rds")

rtracklayer::export.bed(s1.gr,"s1.peaks.bed")

peakAnno_s1f<-annotatePeak(peak = s1.f.gr, TxDb=ens.txdb,
                           tssRegion=c(-3000, 3000),level="gene")

saveRDS(peakAnno_s1f,"peakAnno_s1.f.rds")

rtracklayer::export.bed(s1.f.gr,"s1.peaks.f.bed")
