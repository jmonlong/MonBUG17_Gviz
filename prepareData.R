library(dplyr)
library(magrittr)
library(rtracklayer)
library(GenomicRanges)

## Gencode for chromosome 21
download.file("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz", "gencode.v25.annotation.gtf.gz")
gtf = import("gencode.v25.annotation.gtf.gz")

gencode = gtf %>% as.data.frame %>% filter(seqnames=="chr21", type %in% c("CDS", "UTR", "exon")) %>%
    group_by(transcript_id) %>% filter(all(type!="UTR") | type!="exon") %>%
        mutate(feature=ifelse(type=="UTR", "utr", gene_type)) %>%
            dplyr::rename(symbol=gene_name, transcript=transcript_id) %>%
                select(seqnames, start, end, strand, feature, transcript, symbol) %>%
                    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

save(gencode, file="gencode.v25.ch21.RData")

## Some ranges, e.g. CNVs
cnv = GRanges("chr21", IRanges(runif(10,36400000,37000000), width=runif(10,500,1e4)))
cnv$sample = sample(c("sample1","sample3","sample4"), 10, TRUE)
cnv$type = sample(c("deletion", "duplication"), 10, TRUE)
cnv = c(cnv[1], cnv)
cnv$sample[1] = "sample5"

## QTLs like data
qtls = GRanges("chr21", IRanges(sort(unique(runif(1e3,36400000,37000000))), width=1))
qtls$logPv = -log10(runif(length(qtls)))
qtls$logPv[sample.int(length(qtls),3)] = runif(3, 1,10)

save(cnv, qtls, file="data.RData")

## For GRanges creation example
qtls %>% as.data.frame %>% select(seqnames, start, end, logPv) %>% dplyr::rename(chr=seqnames) %>% write.table(file="qtls.tsv", sep="\t", row.names=FALSE, quote=FALSE)

## Offline backup
library(Gviz)
ideoTrack = IdeogramTrack(genome = "hg19", chromosome = "chr21")
library(AnnotationHub)
ah = AnnotationHub()
query(ah, c("hg19", "dnase"))
dnase = ah[["AH45280"]]
save(ideoTrack, dnase, file="offlineBackup.RData")
