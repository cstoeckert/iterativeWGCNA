#!/usr/bin/env Rscript

## plots eigengene network to a PDF file using WGCNA
## use alterative grDevice to generate other file types (see https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/png.html)

## please view the documentation for the plotEigengeneNetwork
## function (WGCNA) to customize as needed and adjust formatting

require(WGCNA)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    print("Please provide full path to eigengene file name")
    stop("USAGE: eigengene_network_wgcna.R <file>")
}

eigengeneFile <- args[1]
eigengenes <- read.table(eigengeneFile, sep='\t', row.names=1, header=T)

pdf("eigengene-network.pdf", 11, 11)
plotEigengeneNetworks(as.data.frame(t(eigengenes)), c(''), plotDendrograms = TRUE, colorLabels = FALSE, signed = TRUE, plotAdjacency = TRUE)
dev.off()

