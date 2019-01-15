#! /usr/bin/env Rscript
# (c) Konstantin Riege
library("DEXSeq")
library("ggplot2")
library("gplots")
library("BiocParallel")

# samples <- c("/path/to/a1.c","/path/to/a2.c","/path/to/a3.c","/path/to/b1.c","/path/to/b2.c","/path/to/b3.c") 
# conditions <- c("a","a","a","b","b","b")
# labels <- c("a1","a2","a3","b1","b2","b3") 
# libTypes <- c( "paired-end", "paired-end", "paired-end","single-end", "single-end", "single-end")
# gtf <- c("/path/to/dexseqFlattened.gtf")
# threads <- c("48")
# out <- c("/path/to/outdir")
args <- commandArgs(TRUE)

threads <- as.numeric(args[1])
out <- args[2]
samples <- unlist(strsplit(args[3],","))
labels <- unlist(strsplit(args[4],","))
conditions <- unlist(strsplit(args[5],","))
libTypes <- unlist(strsplit(args[6],","))
gff <- args[7]


BPPARAM <- MulticoreParam(workers = threads)

sampleTable <- data.frame(row.names = labels, condition = conditions, libType = libTypes)
ddxDEXSeq <- DEXSeqDataSetFromHTSeq(countfiles = samples, sampleData = sampleTable, design = ~ sample + exon + condition:exon, flattenedfile = gff)

# DEXSeq does:
# ddx <- estimateSizeFactors(ddxDEXSeq)
# ddx <- estimateDispersions(ddx, BPPARAM = BPPARAM)
# plotDispEsts(ddx)
# ddx <- testForDEU(ddx, BPPARAM = BPPARAM)
# ddx <- estimateExonFoldChanges(ddx, BPPARAM = BPPARAM)
# ddxr <- DEXSeqResults(ddx)
ddxr <- DEXSeq(ddxDEXSeq, BPPARAM = BPPARAM)
file <- paste(out, "ddxr.Rdata", sep = "/")
save(ddxr, file = file)

ddxr <- ddxr[order(ddxr$padj) , ]
ddxrNoNA <- ddxr[ !is.na(ddxr$padj) , ]
ddxrNoNA <- ddxrNoNA[ ddxrNoNA$exonBaseMean > 0.0 , ]
ddxrNoNA05 <- ddxrNoNA[ ddxrNoNA$padj < 0.05 , ]
ddxrNoNA01 <- ddxrNoNA[ ddxrNoNA$padj < 0.01 , ]

csv <- paste(out, "full.csv", sep = "/")
o <- as.data.frame(ddxr)
o$transcripts <- ''
# unless removed, output file format gets corrupted
# "ENSG00000000001:E001","ENSG00000000001",...,c("ENST00000000001" "ENST00000000002")
# gets
# "ENSG00000000001:E001","ENSG00000000001",...,c("ENST00000000001" 
# "ENST00000000002")
write.csv(o, file = csv)

csv <- paste(out, "filtered_padjNA.csv", sep = "/")
o <- as.data.frame(ddxrNoNA)
o$transcripts <- ''
write.csv(o, file = csv)

csv <- paste(out, "filtered_padj05.csv", sep = "/")
o <- as.data.frame(ddxrNoNA05)
o$transcripts <- ''
write.csv(o, file = csv)

csv <- paste(out, "filtered_padj01.csv", sep = "/")
o <- as.data.frame(ddxrNoNA01)
o$transcripts <- ''
write.csv(o, file = csv)

pdf(paste(out, "ma_plot.pdf", sep = "/"))
plotMA(ddxr)
graphics.off()

IDs <- row.names(ddxrNoNA05)
uniques <- numeric()
for (i in IDs) {
	uniques <- unique(c(uniques, (strsplit(i, ":")[[1]])[1]))
	if(length(uniques)==50){
		break
	}
}
#uniques <- unique(uniques)[1:min(length(uniques),50)]

out <- file.path(out, "genes")
dir.create(out, showWarnings = FALSE)
for (i in uniques) {
	postscript(file.path(out, paste(i, ".ps", sep = "")))
	plotDEXSeq(ddxr, i, legend = TRUE, lwd = 1)
	graphics.off()
	#postscript(file.path(out, paste(i, "_transcripts.ps", sep = "")))
	#plotDEXSeq(ddxr, i, displayTranscripts = TRUE, legend = TRUE, lwd = 1)
	#graphics.off()
}
