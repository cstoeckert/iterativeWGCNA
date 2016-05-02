#!/usr/bin/env python2.7
"""
uses rpy2 python library to create a namespace for R functions underlying iterativeWGCNA
"""

util_functions = """

# convert numeric data frame to real
numeric2real <- function(df) {
     df * 1.0
}

# wrapper for save object b/c doesn't seem to work with rpy2
saveObject <- function(obj, objName, file) {
   assign(objName, obj)
   save(list=c(objName), file = file)
}

# given WGCNA blocks, extracts and transposes eigengene matrix
# labels columns (samples)
# cleans up module names (removes the "ME")

eigengenes <- function(iteration, blocks, sampleNames) {
    eigengenes <- as.data.frame(t(blocks$MEs))
    colnames(eigengenes) <- sampleNames
    eigengenes <- eigengenes[row.names(eigengenes) != "ME0" & row.names(eigengenes) != "MEgrey", ]
    row.names(eigengenes) <- gsub("ME", paste(iteration, "-", sep=""), row.names(eigengenes))

    eigengenes
}

# given WGCNA blocks and gene names, returns
# a data frame with modules mapped to gene names
modules <- function(blocks, geneNames) {
    as.data.frame(blocks$colors, row.names = geneNames)
}

# extract module members
# does not assume same ordering
extractMembers <- function(module, expr, membership) {
    members = t(membership)[, 1] == module
    expr[members, ]
}


"""

other = """

# update eigengene names and output eigengene file


# process eigensimilarity for each gene, get revised membership
# generate gene lists
# return whether converged (no dropped genes)
evaluateFit <- function(data, blocks, runId, targetDir) {
    d <- as.matrix(data)

    #rpy2 is not acting nicely, so just have to redo this for now
    eigengenes <- t(blocks$MEs)
    row.names(eigengenes) <- gsub("ME", "", row.names(eigengenes))
    colors <- blocks$colors

    esim <- NULL
    rowNames <- NULL

    modules <- NULL

    dropped <- NULL
    kept <- NULL

    for (i in 1:nrow(data)) {
        profile <- d[i, ]
        gene <- row.names(d)[i]
        module <- colors[i]

        if (module != "grey") {
            me <- eigengenes[module, ]
            esim <- c(esim, cor(profile, me))
            rowNames <- c(rowNames, gene)
            modules <- c(modules, paste(runId, module, sep="_"))
            kept <- rbind(kept, data[gene, ])
        }
        else {
            dropped <- rbind(dropped, data[gene, ])
        }

    }

    colnames(kept) <- colnames(data)

    fit <- data.frame(Module = modules, Fit = esim, row.names = rowNames)

    write.table(kept, paste(targetDir, "kept-gene-expression.txt", sep="/"), sep="\t", quote=F)
    write.table(fit, paste(targetDir, "module-membership.txt", sep="/"), sep="\t", quote=F)

    if (!is.null(dropped)) {
        colnames(dropped) <- colnames(data)
        write.table(dropped, paste(targetDir, "dropped-gene-expression.txt", sep="/"), sep="\t", quote=F)
        returnval <- FALSE
    }
    else {
        returnval <- TRUE
    }

    returnval
}

"""


__author__ = "Emily Greenfest-Allen"
__copyright__ = "Copyright 2015, University of Pennsylvania"
