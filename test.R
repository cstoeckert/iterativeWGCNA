suppressPackageStartupMessages(require(WGCNA))
suppressPackageStartupMessages(require(igraph))

# wrapper for save object b/c doesn't seem to work with rpy2
saveObject <- function(obj, objName, file) {
   assign(objName, obj)
   save(list=c(objName), file = file)
}

# calculate an adjacency matrix (uses WGCNA adjacency function)
# and outputs to a file
# depends on WGCNA
calcAdjacencyMatrix <- function(data,
                                 power = 6,
                                 workingDirectory = getwd(),
                                 writeMatrix = FALSE,
                                 filename="fit-adjacency-matrix.txt",
                                 verbose = TRUE
                                 ) {
    if (verbose) {
        print("Calculating Signed Adjacency Matrix for Retained Gene Set (Good Fit)")
    }

    adjacencyMatrix <- adjacency(t(data), power = power, type="signed")

    if (writeMatrix) {
        filename = paste(workingDirectory, filename, sep="/")
        write.table(adjacencyMatrix, filename, sep="\t", quote=FALSE)
    }
    adjacencyMatrix
}

# generates an edge list from an adjacency matrix and
# outputs to file
# depends on igraph write.graph

adj2edgeList <- function(adjMatrix, workingDirectory = getwd(), filename="fit-network-edge-list.txt", verbose=TRUE) {
    if (verbose)
        print("Generating and Exporting Network From Fit Genes")

    g <- graph.adjacency(adjMatrix, mode="undirected", weighted = TRUE)

    filename = paste(workingDirectory, filename, sep="/")
    write.graph(g, filename, format="ncol")

    g
}

# run the blockwiseWGCNA
# see WGCNA manual for parameter details
# depends on WGCNA
bWGCNA <- function(data,
                   outputDir,
                   randomSeed = 12345,
                   minCoreKME = 0.80,
                   power = 6,
                   minModuleSize = 20,
                   minCoreKMESize = 15,
                   maxBlockSize = 5000,
                   networkType = signed,
                   TOMType = signed,
                   deepSplit = 0,
                   allowWGCNAThreads = TRUE,
                   verbose = TRUE) {

    if (verbose)
        wgcnaVerbose = 10
    else
        wgcnaVerbose = 0

    if (allowWGCNAThreads)
        allowWGCNAThreads()

    filename = paste(outputDir, 'blockwiseWGCNA-output.pdf', sep='/')
    pdf(filename)
    # run blockwiseWGCNA
    # data * 1.0 to convert to real (handle the following bug: https://support.bioconductor.org/p/76829/)

    blocks <- blockwiseModules(t(data) * 1.0,
                               randomSeed = randomSeed,
                               power = power,
                               networkType = networkType,
                               minModuleSize = minModuleSize,
                               TOMType = TOMType,
                               minCoreKME = minCoreKME,
                               minCoreKMESize = minCoreKMESize,
                               minKMEtoStay = minCoreKME - 0.05,
                               deepSplit = deepSplit,
                               verbose = wgcnaVerbose)
    dev.off()
    blocks
}

# update eigengene names and output eigengene file
processEigengenes <- function(blocks, sampleNames, runId, targetDir) {
    eigengenes <- t(blocks$MEs)

    colnames(eigengenes) <- sampleNames
    row.names(eigengenes) <- gsub("ME", paste(runId, "_", sep=""), row.names(eigengenes))

    write.table(eigengenes, paste(targetDir, 'eigengenes.txt', sep="/"), quote=F)

    eigengenes
}

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

getDroppedGeneExpression <- function(dataFile, keptDataFile, targetDir) {
    all <- read.table(dataFile, row.names = 1, header=T, sep="\t")
    kept = read.table(keptDataFile, row.names=1, header=T, sep="\t")
    dropped <- all[setdiff(row.names(all), row.names(kept)), ]
    write.table(dropped, paste(targetDir, "pass-initial-gene-expression.txt", sep="/"), quote=F)
}
