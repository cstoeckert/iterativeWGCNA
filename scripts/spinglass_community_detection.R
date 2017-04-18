#!/usr/bin/env Rscript

## performs spinglass community detection on eigengene network
## outputs:
   ### module -> meta-module mapping file and
   ### weighted edge-list defining network (can be imported into cytoscape)
   ### saved R workspace
   ### very ugly R plots of the graph, if outputGraphs flag is set to TRUE (see below)

## please view the documentation for the 
## igraph functions to customize as needed
## and adjust output formatting

## assumes: per-sample replicates are adjacent
## note: edge weights are topological overlap

require(igraph)
require(WGCNA)
require(optparse)


## COMMUNITY DETECTION PARAMETERS
## =================================================
sg.spins <- 100 # 25
sg.start.temp <- 1
sg.stop.temp <- 0.01
sg.cool.fact <- 0.99
sg.gamma <- 2 # 1
sg.gamma.minus <- 1

outputGraphs <- TRUE


## HEATMAP PARAMETERS
## =================================================


byapply <- function(x, by, fun, ...)
{
    # Create index list
    if (length(by) == 1)
    {
        nc <- ncol(x)
        split.index <- rep(1:ceiling(nc / by), each = by, length.out = nc)
    } else # 'by' is a vector of groups
    {
        nc <- length(by)
        split.index <- by
    }
    index.list <- split(seq(from = 1, to = nc), split.index)

    # Pass index list to fun using sapply() and return object
    sapply(index.list, function(i)
            {
                do.call(fun, list(x[, i], ...))
            })
}

scale_rows = function(x){
	m = apply(x, 1, mean, na.rm = T)
	s = apply(x, 1, sd, na.rm = T)
	return((x - m) / s)
}

scale_mat = function(mat, scale){
	if(!(scale %in% c("none", "row", "column"))){
		stop("scale argument shoud take values: 'none', 'row' or 'column'")
	}
	mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
	return(mat)
}


## parse command line args
parseCommandLineArgs <- function() {
    options <- list(make_option(c("-e", "--eigengenes"),
                                help="full path to eigengene file name file name",
                                type = "character"),
                    make_option("--averageReplicates",
                                help="average replicates? to reduce replicate-level variation on metanetwork detection",
                                type="logical",
                                action="store_true",
                                default=FALSE),
                       make_option(c("-n", "--numReplicates"),
                                   help="number of replicates per sample; must provide if --averageReplicates option is set",
                                   type="integer",
                                   default=1
                                   ),
                       make_option(c("-s", "--samples"),
                                   help="comma separated list of labels for samples; for use if --averageReplicates is specified",
                                   type="character"
                                   ),
                    make_option(c("-c", "--corThreshold"),
                                   help="minimum correlation value for edge; if not specified, will filter only on p-value",
                                type="double",
                                default=0.0
                                   )
                     ) # end optionList definition

### TO DO: add validation
    parse_args(OptionParser(option_list = options))
}

## returns a list of three eigengene matrices: original, std, average over replicates
parseEigengenes <- function(filename, numReplicates, averageReplicates, samples) {
    eigengenes <- read.table(filename, sep = "\t", row.names = 1, header = T)
    eigengenes.std <- scale_mat(eigengenes, "column") # scale in case non-scaled eigengene option was used for WGCNA
    eigengenes.avg <- eigengenes.std # to avoid null value
    
    if (numReplicates > 1 && averageReplicates ) {
        eigengenes.avg <- as.data.frame(byapply(eigengenes.std, numReplicates, rowMeans)) # average across replicates
        colnames(eigengenes.avg) <- strsplit(samples, ',')[[1]]
        write.table(eigengenes.avg, "eigengenes-averaged-across-replicates.txt", sep="\t", quote=F)
    }

    list(raw = eigengenes, std = eigengenes.std, avg = eigengenes.avg)
}

## generates the network from a TOM matrix
generateGraph <- function(profiles, corThreshold) {
    cor.and.p.matrix <- corAndPvalue(t(profiles))
    cor.and.p.matrix$c[cor.and.p.matrix$p > 0.05] <- 0 # filter out edges lacking statistical support
    cor.and.p.matrix$c[cor.and.p.matrix$c < corThreshold] <- 0 #filter out by threshold

    adjacency.matrix <- adjacency.fromSimilarity(cor.and.p.matrix$c, type="signed")

    tom <- TOMsimilarity(adjacency.matrix, TOMType="signed")
    ## need to clean up a bit more to get a cleaner graph; set topological overlap < 2 to 0
    tom[tom < 0.2] <- 0
    # diss.tom <- 1.0 - tom
    row.names(tom) <- row.names(profiles)
    colnames(tom) <- row.names(profiles)

    graph <- graph.adjacency(as.matrix(tom), weighted=TRUE, diag=FALSE, mode="max")
    print(paste("# VERTICES:", vcount(graph)))
    print(paste("# EDGES:", ecount(graph)))

    print(paste("IS CONNECTED?" , is.connected(graph)))

    if (outputGraphs) {
        png("meta-network-graph.png")
        plot(graph)
        dev.off()
    }
    
    ## identify and remove unconnected nodes
    isolates <- which(degree(graph, mode = "total") == 0) - 1
   
    if (length(isolates > 0)) {
        print("REMOVING UNCONNECTED NODES")
        print(isolates)
        graph <- delete.vertices(graph, isolates)
        if (outputGraphs) {
            png("meta-network-graph-no-isolates.png")
            plot(graph)
            dev.off()
        }
    }
  
   
    write_graph(graph, "meta-network-graph-edges.txt", format = "ncol")
    graph
}


args <- parseCommandLineArgs()
eigengenes <- parseEigengenes(args$eigengenes, args$numReplicates, args$averageReplicates, args$samples)
graph <- generateGraph(eigengenes$avg, args$corThreshold)

# uses random null-model with same degree distribution
communities <- cluster_spinglass(graph, spins = sg.spins, start.temp=sg.start.temp, stop.temp=sg.stop.temp, cool.fact=sg.cool.fact, update.rule="simple", gamma=sg.gamma, implementation="orig", gamma.minus=sg.gamma.minus)

# see igraph manual for all options for exploring/plotting community data object

m <- membership(communities)

modules <- names(m)
m <- as.numeric(m)
names(m) <- modules

write.table(data.frame("Module"=modules, "Metamodule"=m), "meta-network-communities.txt", sep="\t", quote=F, row.names=FALSE)

s <- sizes(communities)
metamodules <- names(s)
write.table(data.frame("Metamodule"=metamodules,"Size"=as.data.frame(s)$Freq), "meta-module-sizes.txt", sep="\t", quote=F, row.names=FALSE)

if (outputGraphs) {
    png("meta-network-graph-communities.png")
    plot(communities, graph)
    dev.off()
}

save.image("spinglass-workspace.RData")
