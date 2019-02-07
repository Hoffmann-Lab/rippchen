#! /usr/bin/env Rscript
# (c) Jeanne Wilbrand
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(reshape2))
suppressMessages(library(data.table))
args = commandArgs(TRUE)
input = args[1]
output = args[2]

table = read.csv(input, header = TRUE, sep = "\t", dec = ".", stringsAsFactors = FALSE, strip.white = T)

# get cluster count
cluster_count = length(unique(table$ID))

# strip gene and cluster ID
foo = table[,2:(ncol(table)-1)]

# log all data columns
flog = log2(foo+1)

# get mean_pw_ch_log
if (nrow(flog)>1){
    res <- as.data.frame(rowMeans(mapply(function(x,y) abs(x-y), flog[, 2:ncol(flog)], flog[, 1:(ncol(flog)-1)])))
} else {
    res <- as.data.frame(mean(mapply(function(x,y) abs(x-y), flog[, 2:ncol(flog)], flog[, 1:(ncol(flog)-1)])))
}


colnames(res)[1] <- "mean_pw_ch_log"

# append gene name and cluster ID
flog2 = cbind(flog,table[,c(1,ncol(table))], res)

molt <- melt(flog2, id.vars = c("ID", "gene", "mean_pw_ch_log"), measure.vars = colnames(foo))

pdf(output)
#width=((as.integer(sqrt(cluster_count)+1))*2), height = ((as.integer(sqrt(cluster_count)+1))*3))
ggplot(molt, aes(x = variable, y = value, group=gene, color = mean_pw_ch_log)) +
    theme_bw() +
    labs(x = "Type", y = "TPM (+1, log2-scaled)", 
         color = "Mean log2-FC", 
         linetype = "") +
    geom_line(alpha=0.3) +
    stat_summary(aes(group = 1, linetype = ''), fun.y = 'median', geom = 'line', size = 1, 
                 show.legend = TRUE, colour = 'green') +
    scale_color_gradient(low = "blue", high = "red") +
    scale_linetype_discrete(name = "Median of TPMs\nper type") +
    facet_wrap( ~ ID, ncol = (as.integer(sqrt(cluster_count)+1)))
graphics.off()

 
 
 
 
  
