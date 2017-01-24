glibrary(Gviz)

##
## Annotation Track
load("data.RData")
cnv

atrack = AnnotationTrack(cnv, name="CNV")
plotTracks(atrack)

atrack = AnnotationTrack(cnv, name="CNV", feature=cnv$type, deletion="indianred", duplication="steelblue")
plotTracks(atrack)

atrack = AnnotationTrack(cnv, name="CNV", group=cnv$sample, showId=TRUE)
plotTracks(atrack)

atrack = AnnotationTrack(cnv, name="CNV", group=cnv$sample, feature=cnv$type, deletion="indianred", duplication="steelblue", showId=TRUE)
plotTracks(atrack)



##
## Adding more tracks

## Axis and chromosome ideogram
gatrack = GenomeAxisTrack()
ideoTrack = IdeogramTrack(genome = "hg19", chromosome = "chr21")

## A data track: QTLs with pvalues
qtls
qtrack = DataTrack(qtls, data="logPv", name="QTLs", type="h")

## Custom gene model
load("gencode.v25.ch21.RData")
gencode
genetrack = GeneRegionTrack(gencode, genome = "hg19", name = "Gene", transcriptAnnotation="symbol")

## DNAse
library(AnnotationHub)
ah = AnnotationHub()
query(ah, c("hg19", "dnase"))
dnase = ah[["AH45280"]] # UW.Pancreas.ChromatinAccessibility.STL003.DNase.DS20753
dnase
dtrack = DataTrack(dnase, data="score", type="gradient", name="DNase")

## Join tracks
from = 36400000
to =   37000000
plotTracks(list(ideoTrack, gatrack, genetrack, dtrack, atrack, qtrack), from=from, to=to)

## Adjust relative sizes of tracks
plotTracks(list(ideoTrack, gatrack, genetrack, dtrack, atrack, qtrack), from=from, to=to, sizes=c(1,1,4,1,4,4))




## save(ideoTrack, gatrack, genetrack, dtrack, atrack, qtrack, file="demo.RData")
