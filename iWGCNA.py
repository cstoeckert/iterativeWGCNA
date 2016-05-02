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

from collections import OrderedDict
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
stats = importr('stats')
utils = SignatureTranslatedAnonymousPackage(rs.util_functions, "utils")

# ========================
# Utilities
# ========================

def warning(prefix, *objs):
    '''
    wrapper for writing to stderr
    '''
    print(prefix, *objs, file=sys.stderr)

# ========================
# I/O and File Management
# ========================

def create_dir(dirName):
    '''
    check if directory exists in the path, if not create
    '''
    try:
        os.stat(dirName)
    except OSError:
        os.mkdir(dirName)

    return dirName

def write_data_frame(df, fileName, rowLabel):
    '''
    write data frame to file; creates new file
    if none exists, otherwise appends new eigengenes
    to existing file
    '''
    try:
        os.stat(fileName)
    except OSError:
        header = (rowLabel,) + tuple(df.colnames)
        with open(fileName, 'w') as f:
            print("\t".join(header), file=f)
    finally:
        df.to_csvfile(fileName, quote=False, sep="\t", col_names=False, append=True)

def read_data():
    '''
    read gene expression data into a data frame
    and convert numeric (integer) data to real
    '''
    exprData = ro.DataFrame.from_csvfile(CML_ARGS.inputFile, sep='\t', header=True, row_names=1)
    return utils.numeric2real(exprData)

# ========================
# Help & Command Line Args
# ========================

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
                        help="full path to input gene expression file; " +
                        "if full path is not provided, " +
                        "assumes the file is in the working directory",
                        required=True)

    parser.add_argument('-o', '--workingDir',
                        help="R working directory; where output will be saved",
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

    parser.add_argument('--saveBlocks',
                        help="save WGNCA blockwise modules for each iteration",
                        action="store_true")

    parser.add_argument('--saveFinalEigengenesOnly',
                        help="only save final eigengenes; if not specified will " +
                        "save eigengenes for each iteration",
                        action="store_true")


    return parser.parse_args()

# ========================
# R Session Management
# ========================

def initialize_r_workspace():
    '''
    initialize the r working environment
    '''
    if CML_ARGS.allowWGCNAThreads:
        wgcna.allowWGCNAThreads()
    base.setwd(CML_ARGS.workingDir)

# ========================
# iWGNA Functions
# ========================

def calculate_kme(expr, eigengene):
    '''
    calculates eigengene connectivity kme
    between an eigengene and expression data set
    '''
    correlation = base.as_data_frame(stats.cor(base.t(expr), base.t(eigengene)))
    return correlation

def run_wgcna(data, iteration):
    '''
    run wgcna
    '''

    # generate the parameters to the blockwiseModules call
    params = {}
    if CML_ARGS.wgcnaParameters is not None:
        params = CML_ARGS.wgcnaParameters

    params['datExpr'] = base.t(data) # have to transpose before passing to WGCNA
    params['saveTOMFileBase'] = iteration + "-TOM"

    return wgcna.blockwiseModules(**params)

def initialize_membership():
    '''
    initialize membership dictionary
    gene list comes from input data row names (DATA.rownames)
    all genes are initially unclassified
    an ordered dictionary is used to keep values in same order
    as input data
    '''

    membership = OrderedDict((gene, "UNCLASSIFIED") for gene in DATA.rownames)
    return membership

def update_membership(iteration, data, blocks, membership):
    '''
    compares new module membership assignments to
    prexisting ones; updates membership list
    '''
    if membership is None:
        membership = initialize_membership()

    modules = utils.modules(blocks, data.rownames)

    # if the gene is in the subset
    # update, otherwise leave as is
    for gene in data.rownames:
        # .rx returns a FloatVector which introduces
        # a .0 to the numeric labels when converted to string
        # which needs to be removed
        # note: R array starts at index 1, python at 0
        module = str(modules.rx(gene, 1)[0]).replace(".0", "")
        if module in ("0", "grey"):
            module = "UNCLASSIFIED"
        else:
            module = iteration + "-" + module

        membership[gene] = module

    return membership

def write_membership(iteration, membership):
    '''
    writes the membership dictionary to file
    '''
    df = ro.DataFrame(membership)
    df.rownames = (iteration)
    write_data_frame(df, "membership.txt", "Iteration")

def write_eigengenes(ematrix):
    write_data_frame(ematrix, "eigengenes.txt", "Module")

def get_member_expression(module, expr, membership):
    '''
    subsets expression data
    returning expression for only members of the specified module
    '''
    return utils.extractMembers(module, expr, ro.DataFrame(membership))

def set_iteration_label(runId, passId):
	label = "p" if passId == 0 else "r-" + str(passId)
	label = label + "_iter_" + str(runId)
	return label

# ========================
# iWGCNA Main
# ========================

def iWGCNA():
    # initialize iteration counters and iteration labels
    runId = 0
    passId = 0
    iteration = set_iteration_label(runId, passId)

    membership = None
    algConverged = False

    iterationIndex = -1
    data = utils.numeric2real(DATA) #get_expression_subset(DATA, membership, iterationIndex, False)
    warning(data.nrow)

    blocks = run_wgcna(data, iteration)
    if CML_ARGS.saveBlocks:
        utils.saveObject(blocks, "blocks", iteration + "-blocks.RData")

    eigengenes = utils.eigengenes(iteration, blocks, data.colnames)
    write_eigengenes(eigengenes)

    membership = update_membership(iteration, data, blocks, membership)
    write_membership(iteration, membership)

    # eigengene connectivity (kME) calculation
    for module in eigengenes.rownames:
        eg = eigengenes.rx(module, True)
        warning(module)
        moduleMembers = get_member_expression(module, data, membership)
        # warning(membership)
        warning(moduleMembers.nrow)
        warning(moduleMembers)
        # kME = calculate_kme(data, eg)
        # warning(kME)





    # pass convergence: nResiduals = 0
    # alg convergence: nFit = 0

    # restructure: create 2-d matrix of genes by classification per run


    # while not algConverged:
    #     passId = passId + 1

    #     targetDir = args.workingDir + "/pass" + str(passId)
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
        CML_ARGS = parse_command_line_args()
        create_dir(CML_ARGS.workingDir)
        initialize_r_workspace()

        if CML_ARGS.verbose:
            warning("Working directory:", CML_ARGS.workingDir)
            warning("Allowing WGCNA Threads?", "TRUE" if CML_ARGS.allowWGCNAThreads else "FALSE")
            warning("Running WGCNA with the following parameters")
            warning(CML_ARGS.wgcnaParameters)

        DATA = read_data()

        iWGCNA()

    finally:
        warning("iWGCNA: complete")


__author__ = "Emily Greenfest-Allen"
__copyright__ = "Copyright 2015, University of Pennsylvania"
