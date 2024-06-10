library(genomation)
library(readxl)
library(tidyverse)
library(GenomicRanges)

reps<-readGeneric("../annotation/UCSC_mm10_repeat_masker.tsv",chr = 1,start = 2,end = 3,strand=4,keep.all.metadata = T,header = T)

rep_names=c("L1Md_T","L1Md_Gf","L1Md_F","L1Md_A","IAPEy-int","IAPEz-int","MMERVK10C-int")

repList<-rep_names %>% map(function(n){reps[grep(n,reps$repName),]}) %>% set_names(rep_names)
saveRDS(repList,"repList_unfiltered.rds")

## Save rep as bed
library(rtracklayer)
dir.create("bed_files")
names(repList) %>% map(~ export.bed(object=repList[[.]],con=paste0("bed_files/",.,".all.bed")))

### repList filtered by length based on previous analysis
repList[[1]]<-subset(repList[[1]],width(repList[[1]])>=5000)
repList[[2]]<-subset(repList[[2]],width(repList[[2]])>=5000)
repList[[3]]<-subset(repList[[3]],width(repList[[3]])>=5000)
repList[[4]]<-subset(repList[[4]],width(repList[[4]])>=5000)
repList[[5]]<-subset(repList[[5]],width(repList[[5]])>=5000)
repList[[6]]<-subset(repList[[6]],width(repList[[6]])>=5000)
repList[[7]]<-subset(repList[[7]],width(repList[[7]])>=5000)

saveRDS(repList,"repList.rds")
names(repList) %>% map(~ export.bed(object=repList[[.]],con=paste0("bed_files/",.,".filt.bed")))

