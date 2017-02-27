#!/usr/bin/env Rscript

## performs spinglass community detection on eigengene network
## outputs as cytoscape js file

## please view the documentation for the 
## function (igraph) to customize as needed and adjust formatting

require(igraph)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    print("Please provide full path to eigengene file name")
    print("USAGE: eigengene_network_community_detection.R <file>")
    stopifnot(FALSE)
}
