#!/usr/bin/env python2.7
"""
"""

from __future__ import print_function
from __future__ import with_statement

import argparse
import os
import sys


def findBestFit(profileFileName, eigengeneFileName, outputDir):

    annotatedProfiles = pd.read_table(profileFileName, header = 0, index_col = 0)
    annotatedEigengenes = pd.read_table(eigengeneFileName, header = 0, index_col = 0)

    # assume profile file is as follows:
    # gene module samples
    ncols = len(annotatedProfiles.columns.values)

    samples = annotatedProfiles.columns.values[range(1,ncols)]
    modules = annotatedEigengenes.index.values
    genes = annotatedProfiles.index.values

    originalModules = annotatedProfiles["Module"]

    warning("Calculating Fit")
    fit = pd.DataFrame(0, index=genes, columns=modules)

    i = 0
    for g in genes:
        i = i + 1
        if (i % 250 == 0):
            warning("Processing", i, "th gene:", g)

        profile = annotatedProfiles.loc[g, samples].tolist()
        membership = annotatedProfiles.loc[g, "Module"]

        for m in modules:
            eigengene = annotatedEigengenes.loc[m, samples].tolist()
            fit.loc[g,m] = pearsonr(profile, eigengene)[0]

    fit.to_csv(outputDir + "/kme-fit.txt", sep="\t", float_format="%.3f")
    # free up some memory
    annotatedProfiles = []
    annotatedEigengenes = []

    warning("Calculating Fit")

    bestFitModules = fit.idxmax(axis = 1)

    warning("Evaluating Revised Assignments")
    revised = {}
    diff = {}
    for g in genes:
        original = originalModules.loc[g]
        best = bestFitModules.loc[g]

        if (i % 250 == 0):
            warning("Processing", i, "th gene:", g)

        bestFitR = round(fit.loc[g, best], 3)
        if (original == "UNCLASSIFIED"):
            originalR = 0
        else:
            originalR = round(fit.loc[g, original], 3)

        diffR = bestFitR - originalR
        diff[g] = diffR

        if (best == original):
            revised[g] = False
        else:
            if (original == "UNCLASSIFIED"):
                revised[g] = True if (diffR >= 0.8) else False
            elif (originalR < 0.8):
                revised[g] = True if (bestFitR >= 0.8) else False
                bestFitModules.loc[g] = "UNCLASSIFIED" if (bestFitR < 0.8) else best # unclassify things not meeting 0.8 cut-off
            else:
                revised[g] = True if (diffR >= 0.05) else False


    membershipComparison = pd.DataFrame({'BestFitModule' : fit.idxmax(axis = 1), 'OriginalModule' : originalModules, "Difference": diff, "Revised": revised}, index = genes)
    membershipComparison.to_csv(outputDir + "/kme-comparison.txt", sep="\t", float_format="%.3f")

    with open(outputDir + "/kme-shift-module-membership.txt", 'w') as f:
        for g in genes:
            newModule = membershipComparison.loc[g, 'BestFitModule'] if membershipComparison.loc[g, 'Revised'] == True else  membershipComparison.loc[g, 'OriginalModule']
            print(g + '\t' + newModule, file = f)


if __name__ == "__main__":
    cml_parser = argparse.ArgumentParser()
    cml_parser.add_argument('-g', '--geneExpressionFile', help="tab delimited annotated (with module membership as second column) gene expression file with one row per gene and one column per sample (input to itWGCNA)", required=True)
    cml_parser.add_argument('-e', '--eigengeneFile', help="tab delimited file of eigengenes", required=True)
    cml_parser.add_argument('-o', '--outputDir', help="output file dir", required=True)

    args = cml_parser.parse_args()

    findBestFit(args.geneExpressionFile, args.eigengeneFile, args.outputDir)

