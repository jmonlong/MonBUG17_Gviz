library(Gviz)
library(GenomicRanges)

gr = GRanges("chr21", IRanges(runif(3,36400000,37000000), width=c(1e4,10e4,14e4)), name=c("Mon","BUG","Jan"))
atrack = AnnotationTrack(gr, name="Track")
gatrack = GenomeAxisTrack()
qtls = GRanges("chr21", IRanges(sort(unique(runif(1e3,36400000,37000000))), width=1))
qtls$logPv = -log10(runif(length(qtls)))
qtls$logPv[sample.int(length(qtls),3)] = runif(3, 1,10)
qtrack = DataTrack(qtls, name="QTLs")
ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = "chr21")
load("gencode.v25.ch21.RData")
from = 36400000
to =   37000000
genetrack <- GeneRegionTrack(gencode, genome = "hg19", name = "Gene", transcriptAnnotation="symbol")
plotTracks(list(ideoTrack, gatrack, genetrack, atrack, qtrack), from=from, to=to)

save(ideoTrack, gatrack, genetrack, atrack, qtrack, file="demo.RData")
