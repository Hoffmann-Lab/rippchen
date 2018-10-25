#! /usr/bin/env Rscript
# (c) Arne Sahm
library('WGCNA')
library('DGCA')
library('ggplot2')

args = commandArgs(TRUE)
threads = as.numeric(args[1])
memory = as.numeric(args[2])
filter = args[3]
tpms_raw = read.table(args[4], header=T, sep='\t', stringsAsFactors=F, row.names=1)
out = args[5]

#tpms_raw[is.na(tpms_raw)] = 0
#tpms_filtered=tpms_raw[which(rowSums(tpms_raw)!=0),]
#q = c()
#for (i in 1:ncol(tpms_filtered)){
#  q = c(q, unname(quantile(tpms_filtered[,i][tpms_filtered[,i]>0],c(0.05))) )
#}
#tpms=log(tpms_filtered+1)

if (filter) {
    tpms = filterGenes(tpms_raw, filterTypes='central', filterCentralType='median', filterCentralPercentile = 0.3)
} else {
    tpms = tpms_raw
}
tpms = log(tpms+1)
pdf(paste(out, "histogram.pdf", sep = "."))
ggplot(tpms[1], aes(x=tpms[,1])) + 
  geom_histogram(binwidth = 0.05, aes(fill=..count..)) + 
  theme_bw() + 
  theme(legend.box = "horizontal", plot.title = element_text(hjust = 0.5)) + 
  ggtitle('TPM Histogram') + 
  xlab('TPM') + 
  ylab('#Genes')
graphics.off()

enableWGCNAThreads(threads)
sftThreshold = pickSoftThreshold(
  t(tpms),
  dataIsExpr = TRUE,
  weights = NULL,
  RsquaredCut = 0.85,
  powerVector = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
  removeFirst = FALSE,
  nBreaks = 10,
  blockSize = NULL,
  corFnc = bicor,
  corOptions = list(use = 'p'),
  networkType = "signed",
  moreNetworkConcepts = FALSE,
  gcInterval = NULL,
  verbose = 0,
  indent = 0
)

pdf(paste(out, "scale.pdf", sep = "."))
plot(sftThreshold$fitIndices[,1], -sign(sftThreshold$fitIndices[,3])*sftThreshold$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'))
text(sftThreshold$fitIndices[,1], -sign(sftThreshold$fitIndices[,3])*sftThreshold$fitIndices[,2],
     labels=c(seq(1, 10, by = 1), seq(12, 30, by = 2)),cex=1,col='red'); abline(h=0.90,col='red')
graphics.off()

x = -sign(sftThreshold$fitIndices[,3])*sftThreshold$fitIndices[,2]
power = which(x>=0.9)[1]
power = ifelse(!is.na(power),power,ifelse(ncol(tpms)<20,18,ifelse(ncol(tpms)<31,16,ifelse(ncol(tpms)<41,14,12)))) 
wgcna_result = blockwiseModules(
  t(tpms),
  power = power,
  maxPOutliers = 0.1,
  corType = "bicor",
  TOMType = "signed",
  networkType = "signed",
  replaceMissingAdjacencies = T,
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  verbose = 3,
  saveTOMs = T,
  maxBlockSize = blockSize(50000, rectangularBlocks = TRUE, maxMemoryAllocation = 2^31*memory)
)

pdf(paste(out, "modules.pdf", sep = "."))
plotDendroAndColors(main="Gene dendrogram and modules",wgcna_result$dendrograms[[1]], labels2colors(wgcna_result$unmergedColors)[wgcna_result$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
graphics.off()

save(wgcna_result, file=paste(out, "wgcna.Rdata", sep = "."))

write.table(as.data.frame(wgcna_result$colors, row.names=row.names(tpms)), file=paste(out, 'cluster.tsv', sep = "."), quote=FALSE, sep='\t', col.names = F)
write.table(as.data.frame(wgcna_result$unmergedColors, row.names=row.names(tpms)), file=paste(out, 'module.tsv', sep = "."), quote=FALSE, sep='\t', col.names = F)
