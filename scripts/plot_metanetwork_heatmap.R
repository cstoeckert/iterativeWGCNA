#!/usr/bin/env Rscript

## clusters modules within each metamodule
## uses resulting ordering to draw a heatmap

## please view the documentation for the 
## function pheatmap to customize as needed
## to adjust output formatting

## NOTE: pheatmap assigns a limited # of colors to
## annotations (squares indicating metamodule membership
## above each column); thus with > 5-6 or metamodules the
## colors will be cycled.  Use the parameter annotation_colors
## (see ?pheatmap) to provide alternatives

## expects two inputs: eigengene file (e.g., final-eigengenes.txt)
## and membership : tab-delim file with three columns labeled: Gene Module Metamodule

require(pheatmap)
require(optparse)

## SUPPORTING FUNCTIONS
## =================================================

scale_rows <- function(x){
	m <- apply(x, 1, mean, na.rm = T)
	s <- apply(x, 1, sd, na.rm = T)
	(x - m) / s
}

scale_mat <- function(mat, scale){
	if(!(scale %in% c("none", "row", "column"))){
		stop("scale argument shoud take values: 'none', 'row' or 'column'")
	}
	mat <- switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
	mat
}

## parse command line args
parseCommandLineArgs <- function() {
    options <- list(make_option(c("-e", "--eigengenes"), help="full path to eigengene file name file name",
                                   action="store",
                                   type="character"),
                    make_option(c("-m", "--metamodules"), help="full path to module <-> metamodule map",
                                action="store",
                                type="character")
                    ) # end optionList definition
    
    parse_args(OptionParser(option_list = options))
}


## MAIN
## =================================================

args <- parseCommandLineArgs()
eigengenes <- read.table(args$eigengenes, sep = "\t", row.names = 1, header = T)
eigengenes.std <- scale_mat(eigengenes, "column")

membership <- read.table(args$metamodules, sep = "\t", row.names=1, header = T)
metamodules <- unique(membership$Metamodule)

mmProfiles <- NULL
orderedModules <- NULL
## extract modules in each metamodule
## calculate "average" eigengene for the metamodule
for (mm in metamodules) {
    modules <- row.names(membership)[which(membership$Metamodule == mm)]
    mEigengenes <- eigengenes[modules, ]
    mmProfiles <- rbind(mmProfiles, colMeans(mEigengenes))
}

# cluster metamodules
cMatrix <- cor(t(mmProfiles))
h <- hclust(as.dist(1-cMatrix))
orderedMetamodules <- metamodules[h$order]

## cluster modules within each metamodule by their
## eigengenes (hierarchical, using correlation)
## save order
orderedModules <- NULL
for (mm in orderedMetamodules) {
    modules <- row.names(membership)[which(membership$Metamodule == mm)]
    ## print(paste(mm, length(modules)))

    if (length(modules) > 1) {
        mEigengenes <- eigengenes[modules, ]
        cMatrix <- cor(t(mEigengenes))
        h <- hclust(as.dist(1-cMatrix))
        orderedModules <- c(orderedModules, modules[h$order])
    }
    else {
        orderedModules <- c(orderedModules, modules)
    }
}

orderedMembership <- data.frame(row.names=orderedModules, Metamodule=as.factor(membership[orderedModules, ]))
orderedEigengenes <- eigengenes[orderedModules, ]

c <- colorRampPalette(c("green", "green", "black", "red", "red"))(100)

pheatmap(t(orderedEigengenes), scale="column", cluster_cols=F, cluster_rows=F, cellwidth=10, cellheight=10, border_color="black", color=c, annotation_col=orderedMembership,filename="meta-network-heatmap.pdf")
pheatmap(t(orderedEigengenes), scale="column", cluster_cols=F, cluster_rows=F, cellwidth=10, cellheight=10, border_color="black", color=c,annotation_col=orderedMembership, filename="meta-network-heatmap.png")



