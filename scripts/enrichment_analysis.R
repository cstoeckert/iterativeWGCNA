#!/usr/bin/env Rscript

## takes a gene->cluster membership map (tab-delimted)
## and performs gene set enrichment on each cluster
## currently only works for human or mouse, but can be easily customized
## to use alternative annotation databases

require(org.Hs.eg.db)
require(org.Mm.eg.db)
require(GO.db)
require(clusterProfiler)
require(optparse)

## parse command line args
parseCommandLineArgs <- function() {
    options <- list(make_option(c("-f", "--file"), help="gene membership file",
                                   action="store",
                                   type="character"),
                    make_option(c("-m", "--moduleColumnNumber"), help="column # containing cluster (module) membership",
                                action="store",
                                type="integer"),
                    make_option(c("-g", "--geneColumnNumber"), help="column # containing gene identifiers; default assume 1st column",
                                action="store",
                                type="integer", default=1),
                    make_option(c("-i", "--idType"), help="gene id type: entrez or symbol",
                                action="store", default="symbol", 
                                type="character"),
                    make_option(c("-t", "--taxon"), help="taxon: human or mouse",
                                action="store", default="mouse", 
                                type="character"),
                    ) # end optionList definition
    
    parse_args(OptionParser(option_list = options))
}

## MAIN
## =================================================

args <- parseCommandLineArgs()
membership <- read.table(args$file, sep="\t", row.names=args$geneColumnNumber, sep="\t", stringsAsFactors=F)

clusters <- unique(data[, args$moduleColumnNumber])
clusters <- factor(clusters)
data[,args$moduleColumnNumber] <- factor(data[, args$moduleColumnNumber])



data$Entrez <- as.character(data$Entrez)

universe <- data[data$Metamodule != "UNCLASSIFIED", "Entrez"]
universe <- unique(universe)
universe <- universe[!is.na(universe)]


ontologies <- c("BP", "CC", "MF")

for (ontology in ontologies) {
    print(ontology)
    for (m in metaModules) {
        if (m == "UNCLASSIFIED") {
            next
        }

        print(m)
        genes <- data[data$Metamodule == m, "Entrez"]
        genes <- unique(genes[!is.na(genes)])

        level2 <-  groupGO(genes, organism = "mouse", ont = ontology, level = 2, readable = FALSE)
        level3 <-  groupGO(genes, organism = "mouse", ont = ontology, level = 3, readable = FALSE)

        e <- enrichGO(genes
                 , organism = "mouse"
                 , ont = ontology
                 , pvalueCutoff = 0.05
                 , pAdjustMethod = "BH"
                 , universe
                 , qvalueCutoff = 0.2
                 , minGSSize = 5
                 , readable = TRUE
                 )


        write.table(e@result, paste(m, "_", ontology, ".txt", sep=""), sep="\t", quote=F)
    }
}


