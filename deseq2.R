#! /usr/bin/env Rscript
# (c) Konstantin Riege
library("DESeq2")
library("ggplot2")
library("gplots")
library("BiocParallel")
library("genefilter")
library("RColorBrewer")

# threads <- 48
# out <- /path/to/outdir
# samples <- /path/to/a1.c,/path/to/a2.c,/path/to/a3.c,/path/to/b1.c,/path/to/b2.c,/path/to/b3.c 
# labels <- a1,a2,a3,b1,b2,b3
# conditions <- a,a,a,b,b,b
# levels <- a,b
# patients <- a,b,c,a,b,c

args <- commandArgs(TRUE)
threads <- as.numeric(args[1])
out <- args[2]
samples <- unlist(strsplit(args[3],","))
labels <- unlist(strsplit(args[4],","))
conditions <- unlist(strsplit(args[5],","))
levels <- unlist(strsplit(args[6],","))
patients <- unlist(strsplit(args[7],","))

BPPARAM <- MulticoreParam(workers = threads)

if (is.na(patients) || length(unique(patients))==length(patients) ) {
	sampleTable <- data.frame(sampleName = samples, fileName = samples, condition = conditions, type = labels)
	ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = "", design= ~ condition)
} else {
	sampleTable <- data.frame(sampleName = samples, fileName = samples, condition = conditions, type = labels, patients = patients)
	ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = "", design= ~ patients + condition)
}

# avoid internal alphabetic ordering
ddsHTSeq$condition <- factor(ddsHTSeq$condition, levels=levels)
ddsHTSeq$type <- factor(ddsHTSeq$type, levels=labels)

########################
# run DESeq, which does:
# dds <- estimateSizeFactors(ddsHTSeq)
# dds <- estimateDispersions(dds)
# dds <- nbinomWaldTest(dds)
dds <- DESeq(ddsHTSeq, parallel = TRUE, BPPARAM = BPPARAM)
file <- paste(out, "dds.Rdata", sep = "/")
save(dds, file = file)

#rld object for PCA plots
rld <- rlog(dds)
#vsd object for heatmaps
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay(rld)
vstMat <- assay(vsd)
par(mfrow=c(1,3))
notAllZero <- (rowSums(counts(dds)) > 0)

# get and sort results
ddsr <- results(dds)
file <- paste(out, "ddsr.Rdata", sep = "/")
save(ddsr, file = file)
ddsr <- ddsr[rev(order(abs(ddsr$log2FoldChange))) , ]
ddsrNoNA <- ddsr[ !is.na(ddsr$log2FoldChange) , ]
ddsrNoNA <- ddsrNoNA[ !is.na(ddsrNoNA$padj) , ]
ddsrNoNA <- ddsrNoNA[ ddsrNoNA$baseMean > 0.0 , ]
ddsrNoNAp05 <- ddsrNoNA[ ddsrNoNA$padj < 0.05 , ]
ddsrNoNAp01 <- ddsrNoNA[ ddsrNoNA$padj < 0.01 , ]

#######################
# write statistics
normCounts <- counts(dds, normalized = T)
csv <- paste(out, "normalized_counts.csv", sep = "/")
write.csv(as.data.frame(normCounts), file = csv)

csv <- paste(out, "sizeFactors.csv", sep = "/")
write.csv(as.data.frame(dds$sizeFactor), file = csv)

summary(ddsr)
summary <- paste(out, "summary.txt", sep = "/")
cat("#ddsr$padj < 0.1:\nFALSE\tTRUE\n", file = summary)
cat(table(ddsr$padj < 0.1), file = summary, append = TRUE)
cat("\n\n", file = summary, append = TRUE)
cat("#ddsr$padj < 0.05:\nFALSE\tTRUE\n", file = summary, append = TRUE)
cat(table(ddsr$padj < 0.05), file = summary, append = TRUE)
cat("\n", file = summary, append = TRUE)

pdf(paste(out, "dispersion.pdf", sep = "/"))
plotDispEsts(dds)
graphics.off()

pdf(paste(out, "pvalue_histogram.pdf", sep = "/"))
use <- ddsr$baseMean > attr(ddsr, "filterThreshold")
table(use)
h1 <- hist(ddsr$pvalue[!use], breaks = 0:50/50, plot = FALSE)
h2 <- hist(ddsr$pvalue[use], breaks = 0:50/50, plot = FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE, col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)), adj = c(0.5,1.7), xpd=NA)
graphics.off()

#######################
# write tables
csv <- paste(out, "normalized_counts_vsd.csv" ,sep = "/")
write.csv(as.data.frame(assay(vsd)), file = csv)

csv <- paste(out, "full.csv", sep = "/")
write.csv(as.data.frame(ddsr), file = csv)

csv <- paste(out, "filtered_padjNA_fcNA.csv", sep = "/")
write.csv(as.data.frame(ddsrNoNA), file = csv)

csv <- paste(out, "filtered_padj05.csv", sep = "/")
write.csv(as.data.frame(ddsrNoNAp05), file = csv)

csv <- paste(out, "filtered_padj01.csv", sep = "/")
write.csv(as.data.frame(ddsrNoNAp01), file = csv)

#######################
# plots

# MA plot
# These plots show the log2 fold changes from the treatment over the
# mean of normalized counts, i.e. the average of counts normalized by
# size factors. The left plot shows the “unshrunken” log2 fold changes,
# while the right plot, produced by the code above, shows the shrinkage
# of log2 fold changes resulting from the incorporation of zero-centered
# normal prior. The shrinkage is greater for the log2 fold change
# estimates from genes with low counts and high dispersion, as can be
# seen by the narrowing of spread of leftmost points in the right plot.
pdf(paste(out, "ma_plot.pdf", sep = "/"))
plotMA(ddsr)
graphics.off()

# PCA plot 
# to check for batch effects and the like
pdf(paste(out, "pca_plot.pdf", sep = "/"))
data <- plotPCA(rld, intgroup = c("condition", "type"), returnData = TRUE) 
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color = condition, group = condition, shape = labels)) +
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
rownames(data) <- data$sampleNO
ggplot(data, aes(PC1, PC2, color = condition, group = condition, shape = labels)) +
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



# Heatmap
# of the count matrix, top 50 normalized basemean counts
postscript(paste(out, "heatmap_counts_vsd.ps", sep = "/")) # to replace IDs by Names afterwards
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) # "RdBu" for red/blue
topFrame <- as.data.frame(rowMeans(counts(dds, normalized = TRUE)))
colnames(topFrame) <- c("means")
topFrame$lines <- c(1:nrow(topFrame))
topFrame <- topFrame[order(topFrame$means, decreasing=TRUE) , ]
topIDs <- row.names(topFrame)[1:50]  
topLines <- topFrame$lines[1:50]
# Colv = as.dendrogram(hclust(col.dist, method = "centroid"))
heatmap.2(assay(vsd)[topLines , ], col = hmcol, Rowv =TRUE, Colv = FALSE, scale = "none", dendrogram = "row", trace = "none", margin = c(8, 8), labCol = labels, labRow = topIDs, cexRow = 0.6, cexCol = 1, key.title = NA, key.ylab = NA, key.xlab = "vsdCounts")
graphics.off()

# Heatmap
# of the count matrix, top 50 foldchanges
postscript(paste(out, "heatmap_fc_vsd.ps", sep = "/")) # to replace IDs by Names afterwards
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) # "RdBu" for red/blue
IDs <- row.names(dds) # unsorted all from dds
topIDs <- row.names(ddsrNoNA)[1:50] # sorted, no NA, basemean > 0
topLines <- which(IDs %in% topIDs) # grep line numbers of top ids from RAW unsorted all to select from dds
topIDs <- IDs[topLines]
# optional sorting by normalized basemean counts, in case plotting without clustering
topFrame <- as.data.frame(rowMeans(counts(dds, normalized = TRUE))[topLines])
colnames(topFrame) <- c("means")
topFrame$lines <- topLines
topFrame <- topFrame[order(topFrame$means, decreasing=TRUE) , ]
topIDs <- row.names(topFrame)
topLines <- topFrame$lines
# Colv = as.dendrogram(hclust(col.dist, method = "centroid"))
heatmap.2(assay(vsd)[topLines , ], col = hmcol, Rowv =TRUE, Colv = FALSE, scale = "none", dendrogram = "row", trace = "none", margin = c(8, 8), labCol = labels, labRow = topIDs, cexRow = 0.6, cexCol = 1, key.title = NA, key.ylab = NA, key.xlab = "vsdCounts")
graphics.off()


# Heatmap
# of the count matrix, top 50 foldchanges using gene names and GO ids
# library("GOexpress")
# library("biomaRt")
# if (species == "mmu") {
#     ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "dec2015.archive.ensembl.org")
#     biomart <- useDataset("mmusculus_gene_ensembl", ensembl)
# }
# if (species == "hsa") {
#     ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "dec2015.archive.ensembl.org")
#     biomart <- useDataset("hsapiens_gene_ensembl", ensembl)
# }
# IDs <- row.names(ddsr)
# mart <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "name_1006"), filters = "ensembl_gene_id", values = IDs, mart = biomart)
# # GOterm <- "GO:0006955"
# # select <- grep(GOterm, mart$go_id, fixed = TRUE) 
# # GOrld <- rownames(assay(rld)[mart[select , ]$ensembl_gene_id , ])

# pdf(paste(out,"heatmap_fc_vsd.pdf",sep="/"))
# hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) # "RdBu" for red/blue
# IDs <- row.names(ddsr) # unsorted all
# topIDs <- row.names(ddsrNoNA)[1:50] # sorted, no NA, basemean > 0
# topIDsRegex <- numeric()
# for (i in topIDs) {
# 	regex <- (paste("^", i, "$", sep = ""))
# 	topIDsRegex <- c(topIDsRegex, regex)
# }
# regex <- paste(topIDsRegex, collapse = '|')
# topLines <- grepl(regex, IDs) # grep line numbers of top ids from RAW unsorted all to select from dds
# select <- order(topLines, rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:50]
# selectIDs <- row.names(counts(dds, normalized = TRUE)[select , ])
# selectNames <- c()
# for (i in selectIDs) {
# 	select <- grep(i, mart$ensembl_gene_id, fixed = TRUE) 
# 	selectNames <- c(selectNames, rownames(assay(vsd)[mart[select , ]$ensembl_gene_id , ]))
# }
# # Colv = as.dendrogram(hclust(col.dist, method = "centroid"))
# heatmap.2(assay(vsd)[select , ], col = hmcol, Rowv = TRUE, Colv = FALSE, scale = "none", dendrogram = "row", trace = "none", margin = c(8, 6), labCol = labels, labRow = selectNames, cexRow = 0.6, cexCol = 1, key.title = NA, key.ylab = NA, key.xlab = "log2FC")
# graphics.off()
