library(Gviz)
library(GenomicRanges)

## A simple annotation track
gr = GRanges("chr21", IRanges(c(36506580, 36668973, 36906672), width=c(1e4,10e4,14e4)), name=c("Mon","BUG","Jan"))
atrack = AnnotationTrack(gr, name="Track")

## Axis and chromosome ideogram
gatrack = GenomeAxisTrack()
ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = "chr21")

## A data track: QTLs with pvalues
qtls = GRanges("chr21", IRanges(sort(unique(runif(1e3,36400000,37000000))), width=1))
qtls$logPv = -log10(runif(length(qtls)))
qtls$logPv[sample.int(length(qtls),3)] = runif(3, 1,10)
qtrack = DataTrack(qtls, data="logPv", name="QTLs", type="h")

## Custom gene model
load("gencode.v25.ch21.RData")
genetrack <- GeneRegionTrack(gencode, genome = "hg19", name = "Gene", transcriptAnnotation="symbol")

## Join tracks
from = 36400000
to =   37000000
plotTracks(list(ideoTrack, gatrack, genetrack, atrack, qtrack), from=from, to=to)

## Add UCSC track
constrack <- UcscTrack(genome="hg19", chromosome="chr21", track="Conservation", table="phyloP100wayAll", from=from, to=to, trackType="DataTrack", start="start", end="end", data="score", type="hist", window="auto", col.histogram="darkblue", fill.histogram="darkblue", ylim=c(-3.7, 4), name="Conservation")

plotTracks(list(ideoTrack, gatrack, genetrack, constrack, atrack, qtrack), from=from, to=to)

## DNAse
library(AnnotationHub)
ah = AnnotationHub()
query(ah, c("hg19", "dnase"))
dnase = ah[["AH45280"]]
head(dnase)
dtrack = DataTrack(dnase, data="score", type="gradient", name="DNase")

plotTracks(list(ideoTrack, gatrack, genetrack, dtrack, atrack, qtrack), from=from, to=to)

## save(ideoTrack, gatrack, genetrack, atrack, qtrack, dtrack, file="demo.RData")
## load("demo.RData")

library(dplyr)
library(magrittr)

from = 14960001
to = 14985000
load("cnv.RData")
cnv.gr = makeGRangesFromDataFrame(cnv.df, keep.extra.columns=TRUE)
cnv.gr = subsetByOverlaps(cnv.gr, GRanges("chr21", IRanges(from, to)))

cnvtrack = AnnotationTrack(cnv.gr, name="CNV")
plotTracks(cnvtrack, from=from, to=to)

cnvtrack = AnnotationTrack(cnv.gr, name="CNV", group=cnv.gr$sample)
plotTracks(cnvtrack, from=from, to=to, showId=TRUE)

cnvtrack = AnnotationTrack(cnv.gr, name="CNV", group=cnv.gr$sample, feature=cnv.gr$type, deletion="blue", duplication="red")
plotTracks(cnvtrack, from=from, to=to, showId=TRUE)

file = "~/sftp/guillimin/gs/project/mugqic/projects/TWINS_SIMON/GATK_Pipeline.LP6005057-DNA_B02.clean.dedup.recal.bam"
## twin-ch21-deletion-example.bam
bamImport <- function(file, selection){
    require(Rsamtools)
    param = ScanBamParam(flag=scanBamFlag(isProperPair=FALSE), what = c("qname","pos", "qwidth", "strand"), which = selection)
    bam = scanBam(file, param = param)[[1]]
    bam = as.data.frame(bam)
    bam = subset(bam, !is.na(pos))
    bam$chr = seqnames(selection)[1]
    bam$start = bam$pos
    bam$end = bam$start + bam$qwidth
    gr = makeGRangesFromDataFrame(bam, keep.extra.columns=TRUE)
}
alTrack <- AnnotationTrack(gr, group=gr$qname)
plotTracks(alTrack)

alTrack <- AlignmentsTrack("~/sftp/guillimin/gs/project/mugqic/projects/TWINS_SIMON/GATK_Pipeline.LP6005057-DNA_A02.clean.dedup.recal.bam", isPaired=TRUE)
plotTracks(alTrack, chromosome = "21", from = 14864919, to = 14866919)

afrom <- 44945200
ato <- 44947200
alTrack <- AlignmentsTrack(system.file(package = "Gviz", "extdata", "snps.bam"), isPaired = TRUE)
plotTracks(alTrack, chromosome = "chr21", from = afrom, to = ato)
