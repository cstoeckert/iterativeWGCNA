#!/usr/bin/env python2.7
"""
uses rpy2 python library to create a namespace for R functions underlying iterativeWGCNA
"""

FUNCTIONS = """

# convert numeric data frame to real
numeric2real <- function(df) {
     df * 1.0
}

# return 1-matrix
dissMatrix <- function(df) {
    1.0 - df
}

# add a constant value
add <- function(df, value) {
    df + value
}

# return power-weighted matrix
powerWeightMatrix <- function(df, power) {
    df^power
}


# set values < thresshold to 0
filterByThreshold <- function(df, threshold) {
    df[df < threshold] <- 0
    df
}


# set value of matrix diagonal
diag <- function(df, value) {
    diag(df) <- value
    df
}

# wrapper for save object b/c doesn't seem to work with rpy2
saveObject <- function(obj, objName, file) {
   assign(objName, obj)
   save(list=c(objName), file = file)
}

saveBlockResult <- function(blocks, profiles, file) {
   assign('blocks', blocks)
   assign('expression', profiles)
   save(list=c('blocks', 'expression'), file = file)
}

# calculate degree summary for module genes
degree <- function(adjMatrix, members, threshold) {
    adjSubset <- adjMatrix[members, members]
    inDegree = sum(adjSubset >= threshold) / 2
    adjSubset <- adjMatrix[members,  !names(adjMatrix) %in% members]
    outDegree <- sum(adjSubset >= threshold)
    list(kIn=inDegree, kOut=outDegree)
}


# find two closest modules given a similarity threshold
findCloseModules <- function(similarityMatrix, cutHeight) {
     returnVal <- NULL
     d <- 1 - similarityMatrix
     comparison <- d[d > 0 & d <= cutHeight]
     modulesFound <- sum(comparison) > 0
print(cutHeight)

     if (modulesFound) {
         # indexes of closest modules
         indexes <- which(d == min(comparison), arr.ind = TRUE)
         returnVal <- list(m1 = row.names(d)[indexes[1,1]], m2 = row.names(d)[indexes[1,2]], dissimilarity = d[indexes[1,1], indexes[1,2]])
     }
    returnVal
}


# given WGCNA blocks, extracts and transposes eigengene matrix
# labels columns (samples)
# cleans up module names (removes the "ME")

extractEigengenes <- function(iteration, blocks, sampleNames) {
    eigengenes <- as.data.frame(t(blocks$MEs))
    colnames(eigengenes) <- sampleNames
    eigengenes <- eigengenes[row.names(eigengenes) != "ME0" & row.names(eigengenes) != "MEgrey", ]
    row.names(eigengenes) <- gsub("ME", paste(iteration, "_M", sep=""), row.names(eigengenes))
    eigengenes
}


# extract eigengens from list object output
# from moduleEigengenes function
# label columns (samples)
# clean up module names (remove the "ME")
extractRecalculatedEigengenes <- function(elist, sampleNames) {
   eigengenes <- as.data.frame(t(elist$eigengenes))
   colnames(eigengenes) <- sampleNames
    row.names(eigengenes) <- gsub("ME", "", row.names(eigengenes))
    eigengenes
}

# given WGCNA blocks and gene names, returns
# a data frame with modules mapped to gene names
extractModules <- function(blocks, geneNames) {
    as.data.frame(blocks$colors, row.names = geneNames)
}

# extract module members
# does not assume same ordering
extractMembers <- function(module, expr, membership) {
    membership <- unlist(membership)
    membership <- membership[row.names(expr)]
    members <- membership == module
    expr[members, ]
}

# remove unclassified from expression set
removeUnclassified <- function(expr, membership) {
    membership <- unlist(membership)
    membership <- membership[row.names(expr)]
    classified = membership != "UNCLASSIFIED"
    expr[classified, ]
}

# blue, white, red color scale
BlWhRed <- function() {
    colorRampPalette(c("blue", "white", "red"))(100)
}

# white, yellow, red color scale
WhYlRed <- function() {
   colorRampPalette(c("white", "yellow", "red"))(100)
}

# create color scale by passing a string vector of colors
colorScale <- function(colors) {
   colorRampPalette(colors)(100)
}

"""

__author__ = "Emily Greenfest-Allen"
__copyright__ = "Copyright 2016, University of Pennsylvania"
