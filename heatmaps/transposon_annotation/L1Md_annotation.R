library(genomation)
library(GenomicRanges)

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(genomation)

##Get TE annotations
repList<-readRDS("repList.rds")
rep_names=names(repList)

## Get TE meth info from paper
t<-read_tsv("Transposon_Fig6C.tsv",col_names = T)
tList<-rep_names %>% map(function(n){t[grep(n,t$Feature),]}) %>% set_names(rep_names)

## Join anno and meth info
fList<-names(repList) %>% map(~ repList[[.x]] %>% as.data.frame %>% 
                         right_join(tList[[.x]],by=c("seqnames"="Chr","start"="Start","end"="End"))
) %>% set_names(rep_names)


##remove manually excluded f elements
exclude<-readBed("exclude.bed")
fList$L1Md_F<-fList$L1Md_F %>% filter(!(seqnames %in% seqnames(exclude) & start %in% start(exclude)))

## Filter by threshold
young_list=fList %>% map(~ .x %>% mutate(mCG.diff=MeanMiwi2_KO-MeanWT) %>% 
                          filter(!is.nan(mCG.diff),Mismatch_vs_Consensus<=38) %>% 
                          select(-Strand) %>% 
                          makeGRangesFromDataFrame(keep.extra.columns = T)
) 

old_list=fList %>% map(~ .x %>% mutate(mCG.diff=MeanMiwi2_KO-MeanWT) %>% 
                              filter(!is.nan(mCG.diff),Mismatch_vs_Consensus>38) %>%
                              select(-Strand) %>% 
                              makeGRangesFromDataFrame(keep.extra.columns = T)
)

dir.create("bed_files")
names(young_list) %>% map(~ export.bed(object=young_list[[.]],con=paste0("bed_files2/",.,".Y.bed")))
names(old_list) %>% map(~ export.bed(object=old_list[[.]],con=paste0("bed_files2/",.,".O.bed")))

## L1Md annotation based on consensus mismatch
l1to<-readBed("bed_files2/L1Md_T.O.bed")
l1ty<-readBed("bed_files2/L1Md_T.Y.bed")

l1ao<-readBed("bed_files2/L1Md_A.O.bed")
l1ay<-readBed("bed_files2/L1Md_A.Y.bed")

## Monomer annotation
mono_t<-readGeneric("T_mm10.bed",keep.all.metadata = T)
mono_t<-subset(mono_t,!mono_t$V5 %in% c(2,6))
mono_a<-readGeneric("A_mm10.bed",keep.all.metadata = T)
mono_a<-subset(mono_a,!mono_a$V5 %in% c(2,6))

## T
l1to_mono <- l1to[findOverlaps(l1to,mono_t) |> queryHits() |> unique(),] 
l1to_nomono <- l1to[-(findOverlaps(l1to,mono_t) |> queryHits() |> unique()),] 
l1ty_plus <- c(l1ty,l1to_mono) |> sort()

rtracklayer::export.bed(l1ty_plus,"bed_files2/L1Md_T.Y.bed")
rtracklayer::export.bed(l1to_nomono,"bed_files2/L1Md_T.O.bed")

## A
l1ao_mono <- l1ao[findOverlaps(l1ao,mono_a) |> queryHits() |> unique(),] 
l1ao_nomono <- l1ao[-(findOverlaps(l1ao,mono_a) |> queryHits() |> unique()),] 
l1ay_plus <- c(l1ay,l1ao_mono) |> sort()

rtracklayer::export.bed(l1ay_plus,"bed_files2/L1Md_A.Y.bed")
rtracklayer::export.bed(l1ao_nomono,"bed_files2/L1Md_A.O.bed")
