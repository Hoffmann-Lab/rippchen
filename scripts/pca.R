#! /usr/bin/env Rscript
# (c) Konstantin Riege
suppressMessages(library("DESeq2"))
suppressMessages(library("ggplot2"))
suppressMessages(library("gplots"))
suppressMessages(library("BiocParallel"))
suppressMessages(library("genefilter"))
suppressMessages(library("RColorBrewer"))

# threads <- 48
# out <- /path/to/outdir
# samples <- /path/to/a1.c,/path/to/a2.c,/path/to/a3.c,/path/to/b1.c,/path/to/b2.c,/path/to/b3.c 
# labels <- a1,a2,a3,b1,b2,b3
# replicates <- 1,2,3,1,2,3
# conditions <- a,a,a,b,b,b
# levels <- a,b
# patients <- a,b,c,a,b,c

args <- commandArgs(TRUE)
threads <- as.numeric(args[1])
out <- args[2]
samples <- unlist(strsplit(args[3],","))
labels <- unlist(strsplit(args[4],","))
replicates <- unlist(strsplit(args[5],","))
conditions <- unlist(strsplit(args[6],","))
levels <- unlist(strsplit(args[7],","))
patients <- unlist(strsplit(args[8],","))

BPPARAM <- MulticoreParam(workers = threads)

sampleTable <- data.frame(sampleName = samples, fileName = samples, condition = conditions, type = labels)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = "", design= ~ condition)

ddsHTSeq$condition <- factor(ddsHTSeq$condition, levels=levels)
ddsHTSeq$type <- factor(ddsHTSeq$type, levels=labels)

dds <- DESeq(ddsHTSeq, parallel = TRUE, BPPARAM = BPPARAM)
file <- paste(out, "dds.Rdata", sep = "/")
save(dds, file = file)

rld <- rlog(dds)

pdf(paste(out, "pca_plot.pdf", sep = "/"))
data <- plotPCA(rld, intgroup = c("condition", "type"), returnData = TRUE)
if (is.na(patients) || length(unique(patients))==length(patients) ) {
	data$species <- labels
} else {
	data$species <- patients
}

percentVar <- round(100 * attr(data, "percentVar"))
# ggplot(data, aes(PC1, PC2, color = condition, shape = species)) +
ggplot(data, aes(PC1, PC2, color = condition, group = condition, shape = replicates)) +
	ggtitle(paste("PC1 vs PC2: ", length(rownames(rld)), " genes")) +
	scale_shape_manual(values=1:length(labels)) +
	coord_fixed() + 
	theme_bw() +
	theme(legend.box = "horizontal") +
	geom_point(size = 3) +
	stat_ellipse() +
	xlab(paste0("PC1: ", percentVar[1], "% variance")) +
	ylab(paste0("PC2: ", percentVar[2], "% variance"))
	# ggsave(paste(out, "pca_plot.svg", sep = "/"))
graphics.off()

pdf(paste(out, "pca_plot_top500.pdf", sep = "/"))
topN <- 500
Pvars <- rowVars(assay(rld))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(topN, length(Pvars)))]
PCA <- prcomp(t(assay(rld)[select, ]), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
data <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], PC3 = PCA$x[,3], PC4 = PCA$x[,4], sampleNO = colData(rld)$type, condition = colData(rld)$condition)
if (is.na(patients) || length(unique(patients))==length(patients) ) {
	data$species <- labels
} else {
	data$species <- patients
}
rownames(data) <- data$sampleNO
# ggplot(data, aes(PC1, PC2, color = condition, shape = species)) +
ggplot(data, aes(PC1, PC2, color = condition, group = condition, shape = replicates)) +
	ggtitle(paste("PC1 vs PC2: top ", topN, " variable genes")) +
	scale_shape_manual(values=1:length(labels)) +
	coord_fixed() + 
	theme_bw() +
	theme(legend.box = "horizontal") +
    geom_point(size = 3) +
    stat_ellipse() +
    xlab(paste0("PC1: ",percentVar[1], "% variance")) +
    ylab(paste0("PC2: ",percentVar[2], "% variance"))
	# ggsave(paste(out, "pca_plot_top500.svg", sep = "/"))
graphics.off()
