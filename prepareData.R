library(dplyr)
library(magrittr)
library(rtracklayer)


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
