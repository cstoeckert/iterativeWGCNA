#!/usr/bin/env python2.7
"""
Perform iterative WGCNA analysis

python dependencies:
  * rpy2

R dependencies:
  * igraph
  * WGCNA
"""

from __future__ import print_function
from __future__ import with_statement

import argparse
import os
import sys
import rpy2.robjects as ro
from rpy2.robjects.packages import importr, SignatureTranslatedAnonymousPackage
import rSnippets as rs

# R Package Imports
# ========================

wgcna = importr('WGCNA')
igraph = importr('igraph')
base = importr('base')
utils = SignatureTranslatedAnonymousPackage(rs.util_functions, "utils")


# Parameter Types
# ========================
def restricted_float(value):
    '''
        for argument parsing; restricts float value from 0 to 1
    '''
    value = float(value)
    if value < 0.0 or value > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(value,))
    return value

def parameter_list(strValue):
    '''
    for argument parsing;
    converts a comma separated list of 'param=value' pairs
    into a parm:value hash
    '''

    params = {}
    pairs = strValue.split(',')
    for p in pairs:
        name, value = p.split('=')

        if value in ['TRUE', 'T', 'True', 't']:
            value = True
        if value in ['FALSE', 'F', 'False', 'f']:
            value = False

        params[name] = value

    # if a preference for module labels is not provided
    # by the user, set WGCNA to return numeric labels
    if "numericLabels" not in params:
        params["numericLabels"] = True

    return params

# Utilities
# ========================
def warning(prefix, *objs):
    '''
    wrapper for writing to stderr
    '''
    print(prefix, *objs, file=sys.stderr)

def create_dir(dirName):
    '''
    check if directory exists in the path, if not create
    '''
    try:
        os.stat(dirName)
    except OSError:
        os.mkdir(dirName)

    return dirName


def helpEpilog():
    '''
    text for help epilog
    '''

    inputFileFormatHelp = """ INPUT FILE FORMAT: TBA """
    wgcnaParametersHelp = """ PARAMETERS: DESCRIBE FORMATS AND DEFAULTS """

    return inputFileFormatHelp + "\n\n" + wgcnaParametersHelp

def parse_command_line_args():
    '''
        parse command line args
        '''

    parser = argparse.ArgumentParser(prog="iWGCNA",
                                     description="perform interative WGCNA analysis",
                                     epilog=helpEpilog(),
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--inputFile',
                        metavar='<gene expression file>',
                        help="full path to input gene expression file",
                        required=True)

    parser.add_argument('-o', '--outputFilePath',
                        help="target directory for output files",
                        metavar='<output dir>',
                        default=os.getcwd())

    parser.add_argument('-v', '--verbose',
                        help='print status messages',
                        action="store_true")

    parser.add_argument('-p', '--wgcnaParameters',
                        help="comma separated list of parameters to be passed to WGCNA's " +
                        "blockwiseModules function; " +
                        "e.g., power=6,randomSeed=1234875; " +
                        " see WGCNA manual",
                        type=parameter_list)

    parser.add_argument('--allowWGCNAThreads',
                        action="store_true")


    return parser.parse_args()

def read_data(inputFile):
    '''
    read gene expression data into a data frame
    and convert numeric (integer) data to real
    '''
    exprData = ro.DataFrame.from_csvfile(inputFile, sep='\t', header=True, row_names=1)
    return utils.numeric2real(exprData)

def initialize_r_workspace(args):
    '''
    initialize the r working environment
    '''
    if args.allowWGCNAThreads:
        wgcna.allowWGCNAThreads()
    base.setwd(args.outputFilePath)

# iWGCNA
# ========================

def save_r_object(data, objName, path, fileName):
    '''
    (rename) and save an R object to the specified file path as an .RData file
    :param data: the R object
    :param objName: the name of the object in the save
    :param path: full path to the file
    :param fileName: name of the save
    '''

    fileName = path + '/' + fileName
    rs.wgcna.saveObject(data, objName, fileName)

def write_eigengenes(blocks, samples, runId, outputDirPath):
    '''
    write eigengenes to file
    '''
    return rs.wgcna.processEigengenes(blocks, samples, runId, outputDirPath)

def process_blocks(data, blocks, runId, targetDir, verbose):
    '''
    evaluate blocks by extracting eigengenes,
    calculating eigengene similarities,
    and determining goodness of fit results
    '''

    if verbose:
        warning("Processing BlockwiseResult")

    if verbose:
        warning("Writing Eigengenes")

    eigengenes = write_eigengenes(blocks, data.colnames, runId, targetDir)

    algConverged = False if eigengenes.nrow > 1 else True # no modules detected

    if verbose and passConverged:
        warning("pass converged")

    if not passConverged:
        if verbose:
            warning("Evaluating Fit")
        runConverged = rs.wgcna.evaluateFit(data, blocks, runId, targetDir)

    return (runConverged, algConverged)

def process_run(data, iteration, args):
    '''
    run wgcna
    '''

    # generate the parameters to the blockwiseModules call
    params = {}
    if args.wgcnaParameters is not None:
        params = args.wgcnaParameters

    params['datExpr'] = base.t(data) # have to transpose before passing to WGCNA
    params['saveTOMFileBase'] = iteration

    blocks = wgcna.blockwiseModules(**params)
    warning(blocks)

def get_expression_subset(expr, result, index, isClassified):
    if result is None:
        return expr
    elif isClassified:
        return expr[result[index] != 'UNCLASSIFIED', ]
    else:
        return expr[result[index] == 'UNLASSIFIED',]

def iWGCNA(args):
    runId = 0
    passId = 0
    algConverged = False

    iterationIndex = -1
    iteration = "initial" if passId == 0 else "residual_" + str(passId)
    iteration = iteration + "_" + str(runId)

    exprData = read_data(args.inputFile)



    result = None

    data = get_expression_subset(exprData, result, iterationIndex, False)
    process_run(data, iteration, args)


    # pass convergence: nResiduals = 0
    # alg convergence: nFit = 0

    # restructure: create 2-d matrix of genes by classification per run


    # while not algConverged:
    #     passId = passId + 1

    #     targetDir = args.outputFilePath + "/pass" + str(passId)
    #     create_dir(targetDir)

    #     # first pass -> return all data; otherwise get residuals
    #     data = get_expression_subset(exprData, result, iterationIndex, False)

    #     passConverged = False
    #     runId = 0
    #     while not passConverged:
    #         iterationIndex = iterationIndex + 1

    #         runId = runId + 1
    #         targetDir = targetDir + "/run" + str(runId)
    #         create_dir(targetDir)
    #         key = "pass" + str(passId) + "_run" + str(runId)

    #         passConverged, algConverged = process_run(data, key, targetDir, args)

    #                     # set data to classified genes
    #         data = get_expression_subset(exprData, result, iterationIndex, True)

if __name__ == "__main__":

    try:
        cml_args = parse_command_line_args()
        create_dir(cml_args.outputFilePath)
        initialize_r_workspace(cml_args)

        if cml_args.verbose:
            warning("Working directory:", cml_args.outputFilePath)
            warning("Allowing WGCNA Threads?", "TRUE" if cml_args.allowWGCNAThreads else "FALSE")
            warning("Running WGCNA with the following parameters")
            warning(cml_args.wgcnaParameters)

        iWGCNA(cml_args)

    finally:
        warning("iWGCNA: complete")


__author__ = "Emily Greenfest-Allen"
__copyright__ = "Copyright 2015, University of Pennsylvania"
