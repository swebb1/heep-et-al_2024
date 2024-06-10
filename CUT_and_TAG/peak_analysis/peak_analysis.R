library(GenomicRanges)
library(genomation)
library(ChIPseeker)
library(stringr)
library(dplyr)
library(readr)
library(purrr)
library(ggplot2)

## Read in combined peaks
k4f <- list.files("../",pattern = "H3K4me3_E145_all_peaks.narrowPeak",full.names = T,recursive = F)
names(k4f) <- k4f |> str_remove("../") |> str_remove("_peaks.narrowPeak")
k4.gr <- readPeakFile(k4f)

k9f <- list.files("../",pattern = "H3K9me3_E145_all_peaks.narrowPeak",full.names = T,recursive = F)
names(k9f) <- k9f |> str_remove("../") |> str_remove("_peaks.narrowPeak")
k9.gr <- readPeakFile(k9f)

## Read in individual peaks
k4f <- list.files("../",pattern = "H3K4me3_E145_R[1-2]_all_peaks.narrowPeak",full.names = T,recursive = F)
names(k4f) <- k4f |> str_remove("../") |> str_remove("_all_peaks.narrowPeak")
k4.grl <- k4f |> map(~readPeakFile(.x))
names(k4.grl) <- names(k4f)

k9f <- list.files("../",pattern = "H3K9me3_E145_R[1-2]_all_peaks.narrowPeak",full.names = T,recursive = F)
names(k9f) <- k9f |> str_remove("../") |> str_remove("_all_peaks.narrowPeak")
k9.grl <- k9f |> map(~readPeakFile(.x))
names(k9.grl) <- names(k4f)

# Subset peaks by score, enrichment and no of reads
pt = 86 # Peak threshold
rct = 20

k4rc <- read_tsv("../H3K4me3_E145_peaks_sum.all.tsv",col_names = c("chr","start","end","H3K4me3_E145_R1","H3K4me3_E145_R2"),skip = 1) |>
  mutate(start=start+1,length=end-start)

k9rc <- read_tsv("../H3K9me3_E145_peaks_sum.all.tsv",col_names = c("chr","start","end","H3K9me3_E145_R1","H3K9me3_E145_R2"),skip = 1) |>
  mutate(start=start+1,length=end-start)

rck4 = left_join(as.data.frame(k4.gr),k4rc,by = join_by(seqnames==chr,start==start,end==end)) |>
  filter(H3K4me3_E145_R1 >= rct) |> pull(V4)

rck9 = left_join(as.data.frame(k9.gr),k9rc,by = join_by(seqnames==chr,start==start,end==end)) |>
  filter(H3K9me3_E145_R1 >= rct) |> pull(V4)

k4.f.gr <- subset(k4.gr,k4.gr$V5 >= pt & k4.gr$V4 %in% rck4 & width(k4.gr) < 8000)
k9.f.gr <- subset(k9.gr,k9.gr$V5 >= pt & k9.gr$V4 %in% rck9)

filt_overlaps <- function(gr,grl,overlap=2){
  gr$hits = 0
  for(i in 1:length(grl)){
    qh = findOverlaps(gr,grl[[i]]) |> queryHits()
    gr$hits[qh] = gr$hits[qh] + 1 
  }
  subset(gr,gr$hits>=overlap)
}

k4.fo.gr <- filt_overlaps(k4.f.gr,k4.grl)
k9.fo.gr <- filt_overlaps(k9.f.gr,k9.grl)

co.gr <- GenomicRanges::intersect(k4.fo.gr,k9.fo.gr)
co.gr <- subset(co.gr, width(co.gr) >= 100) #Minimum of 100bases

co.gr$name <- paste0("co_",1:length(co.gr))

seqlevelsStyle(co.gr) <- "UCSC"
saveRDS(co.gr,"co.peaks.rds")

## Annotate regions
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

txdb<-TxDb.Mmusculus.UCSC.mm10.knownGene
ens.txdb<-makeTxDbFromGFF("Mus_musculus.GRCm38.78.edited.gtf",format = "gtf",organism = "Mus musculus")

peakAnno_co<-annotatePeak(peak = co.gr, TxDb=ens.txdb,
                          tssRegion=c(-3000, 3000),level="gene")

saveRDS(peakAnno_co,"peakAnno_co.rds")

rtracklayer::export.bed(co.gr,"co.peaks.bed")
