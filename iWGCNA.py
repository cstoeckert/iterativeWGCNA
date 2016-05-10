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
import logging
from time import strftime
import rpy2.robjects as ro
from rpy2.robjects.packages import importr, SignatureTranslatedAnonymousPackage
import rSnippets as rs

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
            print('\t'.join(header), file=f)
    finally:
        df.to_csvfile(fileName, quote=False, sep='\t', col_names=False, append=True)

def read_data(fileName):
    '''
    read gene expression data into a data frame
    and convert numeric (integer) data to real
    '''
    try:
        data = ro.DataFrame.from_csvfile(fileName, sep='\t', header=True, row_names=1)
    except:
        logger.error("Unable to open input file: " + fileName)
        sys.exit(1)
    return utils.numeric2real(data)

def transpose_file_contents(fileName, rowLabel):
    '''
    read in a file to a dataframe, transpose, and output
    '''
    contents = read_data(fileName)
    contents = ro.DataFrame.transpose(contents)
    write_data_frame(contents, fileName, rowLabel)

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

    return params

def helpEpilog():
    '''
    text for help epilog
    '''

    inputFileFormatHelp = ''' 
------------------------------------------------------
Input File Format
------------------------------------------------------
iWGCNA expects a tab-delimited text file containing 
gene expression data arranged such that there is one
row per gene and one column per sample.  The first column
should contain unique gene identifiers.  For example:

GENE    Sample1    Sample2    Sample3
Gata1    500    715    1000
Phtf2    60    1000    1600
'''
    wgcnaParametersHelp = '''
------------------------------------------------------
WGCNA Parameters
------------------------------------------------------
iWGCNA can accept any parameter valid for the WGCNA
blockwiseModules function.  

See http://www.inside-r.org/packages/cran/wgcna/docs/blockwiseModules

To specify these parameters use the --wgcnaParameters flag followed by
a comma separated list of parameter=value pairs.  For example:

--wgcnaParameters maxBlockSize = 5000,corType=bicor,power=10

sets the maximum block size to 5000 genes,
the correlation type to the biweight correlation,
and the power-law scaling factor (beta) to 10

If parameters are not specified, iWGCNA uses the default WGCNA settings,
except for the following:

minModuleSize=20
saveTOMs=TRUE
minKMEtoStay=0.8
minCoreKME=0.8
reassignThreshold=0.05
networkType=signed
numericLabels=TRUE

'''

    return inputFileFormatHelp + '\n\n' + wgcnaParametersHelp

def parse_command_line_args():
    '''
    parse command line args
    '''

    parser = argparse.ArgumentParser(prog='iWGCNA',
                                     description="perform interative WGCNA analysis",
                                     epilog=helpEpilog(),
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--inputFile',
                        metavar='<gene expression file>',
                        help="full path to input gene expression file; "
                        + "if full path is not provided,\n"
                        + "assumes the file is in the working directory\n;"
                        + "see below for input file format",
                        required=True)

    parser.add_argument('-o', '--workingDir',
                        help="R working directory; where output will be saved",
                        metavar='<output dir>',
                        default=os.getcwd())

    parser.add_argument('-v', '--verbose',
                        help="print status messages",
                        action='store_true')

    parser.add_argument('-p', '--wgcnaParameters',
                        metavar='<param list>',
                        help="comma separated list of parameters to be passed to WGCNA's "
                        + "blockwiseModules function\n"
                        + "e.g., power=6,randomSeed=1234875\n"
                        + "see WGCNA manual & more info below",
                        type=parameter_list)

    parser.add_argument('--allowWGCNAThreads',
                        help="allow WGCNA to use threading;\nsee WGCNA manual",
                        action='store_true')

    parser.add_argument('--saveBlocks',
                        help="save WGNCA blockwise modules for each iteration",
                        action='store_true')

    args = parser.parse_args()
    args.wgcnaParameters = set_wgcna_parameter_defaults(args.wgcnaParameters)

    return args

def set_wgcna_parameter_defaults(params):
    '''
    set default values for WGCNA blockwiseModules
    numericLabels = TRUE
    minKMEtoStay = 0.8
    '''

    if params is None:
        params = {}

    if 'numericLabels' not in params:
        params['numericLabels'] = True
    if 'minKMEtoStay' not in params:
        params['minKMEtoStay'] = 0.8
    if 'minCoreKME' not in params:
        params['minCoreKME'] = 0.8
    if 'minModuleSize' not in params:
        params['minModuleSize'] = 20
    if 'reassignThreshold' not in params:
        params['reassignThreshold'] = 0.05 # 0.0000001 # 1e-6

    return params


# ========================
# Logging
# ========================

def initialize_log():
    '''
    initialize log by setting path and file format
    '''
    logging.basicConfig(filename=os.path.join(CML_ARGS.workingDir, 'iwgcna.log'),
                        filemode='w', format='%(levelname)s: %(message)s',
                        level=logging.DEBUG)

    return logging.getLogger(__name__)

def log_parameters():
    '''
    log WGCNA parameter choices and working directory name
    '''
    logger.info(strftime("%c"))
    logger.info("Working directory: " + CML_ARGS.workingDir)
    logger.info("Allowing WGCNA Threads? "
                + ("TRUE" if CML_ARGS.allowWGCNAThreads else "FALSE"))
    logger.info("Running WGCNA with the following parameters:")
    logger.info(CML_ARGS.wgcnaParameters)

    if CML_ARGS.verbose:
        warning("Working directory: " + CML_ARGS.workingDir)
        warning("Allowing WGCNA Threads? "
                + ("TRUE" if CML_ARGS.allowWGCNAThreads else "FALSE"))
        warning("Running WGCNA with the following parameters:")
        warning(CML_ARGS.wgcnaParameters)

def log_input_data():
    '''
    log size of input data
    '''
    logger.info("INPUT DATA: " + CML_ARGS.inputFile)
    logger.info(str(DATA.ncol) + " SAMPLES")
    logger.info(str(DATA.nrow) + " GENES")

    if CML_ARGS.verbose:
        warning("INPUT DATA: " + CML_ARGS.inputFile)
        warning(str(DATA.ncol) + " SAMPLES")
        warning(str(DATA.nrow) + " GENES")

def log_pass_completion(passId, iteration, initialGeneCount,
                        classifiedGeneCount, moduleCount):
    '''
    summarize pass in log when convergence condition is met
    '''
    message = "Pass " + str(passId) + " converged on iteration " + iteration + "."
    message = message + " FIT: " + str(classifiedGeneCount) + "; RESIDUAL: "
    message = message + str(initialGeneCount - classifiedGeneCount)
    message = message + "; MODULES: " + str(moduleCount)
    logger.info(message)

# ========================
# R Session Management
# ========================

def initialize_r_workspace():
    '''
    initialize the r working environment
    '''
    # set working directory
    base.setwd(CML_ARGS.workingDir)

    # suppress warnings
    ro.r['options'](warn=-1)

    # allow WGCNA to run mutli-threaded
    if CML_ARGS.allowWGCNAThreads:
        wgcna.allowWGCNAThreads()

# ========================
# iWGNA Functions
# ========================

def calculate_kme(expr, eigengene, calculateP):
    '''
    calculates eigengene connectivity kme
    between an eigengene and expression data set
    '''
    if calculateP:
        correlation = wgcna.corAndPvalue(base.t(expr), base.t(eigengene))
    else:
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
    params['saveTOMFileBase'] = iteration + '-TOM'

    return wgcna.blockwiseModules(**params)

def initialize_kme():
    '''
    initialized eigengene connectivity (kME)
    dictionary
    gene list comes from input data row names (DATA.rownames)
    all kMEs are initially NaN
    an ordered dictionary is used to keep values in the same
    order as input data
    '''
    kME = OrderedDict((gene, float('NaN')) for gene in DATA.rownames)
    return kME

def update_kme(kME, data, membership, eigengenes):
    '''
    finds new module membership and updates eigengene
    connectivity (kME)
    for each module, extracts the member subset from the
    expression data and calculates the kME between the module
    eigengene and each member
    '''

    if kME is None:
        kME = initialize_kme()

    for module in eigengenes.rownames:
        moduleEigengene = eigengenes.rx(module, True)
        moduleMemberExpression = get_member_expression(module, data, membership)
        memberKME = calculate_kme(moduleMemberExpression, moduleEigengene, False)
        for gene in memberKME.rownames:
            kME[gene] = memberKME.rx(gene, 1)[0]

    return kME

def write_kme(iteration, kME):
    '''
    writes the eigengene connectivity (kME)
    dictionary to file
    '''
    df = ro.DataFrame(kME)
    df.rownames = (iteration)
    write_data_frame(df, 'eigengene-connectivity.txt', 'Iteration')

def initialize_membership():
    '''
    initialize membership dictionary
    gene list comes from input data row names (DATA.rownames)
    all genes are initially unclassified
    an ordered dictionary is used to keep values in same order
    as input data
    '''

    membership = OrderedDict((gene, 'UNCLASSIFIED') for gene in DATA.rownames)
    return membership

def update_membership(iteration, genes, blocks, membership):
    '''
    compares new module membership assignments to
    prexisting ones; updates membership list
    '''
    if membership is None:
        membership = initialize_membership()

    modules = utils.modules(blocks, genes)

    # if the gene is in the subset
    # update, otherwise leave as is
    for g in genes:
        # .rx returns a FloatVector which introduces
        # a .0 to the numeric labels when converted to string
        # which needs to be removed
        # note: R array starts at index 1, python at 0
        module = str(modules.rx(g, 1)[0]).replace('.0', '')
        if module in ('0', 'grey'):
            module = 'UNCLASSIFIED'
        else:
            module = iteration + '-' + module

        membership[g] = module

    return membership

def write_membership(iteration, membership, isPruned):
    '''
    writes the membership dictionary to file
    :param iteration    iWGCNA iteratoin
    :param membership   gene->module mapping
    :param initialClassificaton    boolean flag indicating whether pruning has been done
    '''
    df = ro.DataFrame(membership)
    df.rownames = (iteration)
    fileName = 'pruned-membership.txt' if isPruned else 'initial-membership.txt'

    write_data_frame(df, fileName, 'Iteration')

def write_gene_counts(iteration, initial, fit):
    '''
    writes the number of kept and dropped genes at the end of an iteration
    '''
    fileName = 'iteration-gene-count-summary.txt'
    try:
        os.stat(fileName)
    except OSError:
        header = ('Iteration', 'Initial', 'Fit', 'Residual')
        with open(fileName, 'a') as f:
            print('\t'.join(header), file=f)
    finally:
        with open(fileName, 'a') as f:
            print('\t'.join((iteration, str(initial),
                             str(fit), str(initial - fit))), file=f)

def write_eigengenes(ematrix):
    '''
    writes the eigengene matrix to file
    '''
    write_data_frame(ematrix, 'eigengenes.txt', 'Module')

def get_member_expression(module, expr, membership):
    '''
    subsets expression data
    returning expression for only members of the specified module
    '''
    return utils.extractMembers(module, expr, ro.ListVector(membership))

def get_residuals(expr, membership):
    '''
    subsets expression data
    returning expression for only residuals to the fit
    '''
    if membership is None:
        return expr
    else:
        return get_member_expression('UNCLASSIFIED', expr, membership)

def remove_residuals(expr, membership):
    '''
    subsets expression data
    removing residuals to the fit
    '''
    if membership is None:
        return expr
    else:
        return utils.removeUnclassified(expr, ro.ListVector(membership))# ro.DataFrame(membership))

def set_iteration_label(runId, passId):
    '''
    generates the unique label for the iteration
    '''
    label = 'p' if passId == 0 else 'r' + str(passId)
    label = label + '_iter_' + str(runId)
    return label

def remove_small_modules(memberCount, membership, genes):
    '''
    checks membership counts and removes
    any modules that are too small
    by updating gene membership to UNCLASSIFIED
    returns updated membership and count of total number of modules
    as well as count of total classified genes
    '''
    modules = {}
    classifiedGeneCount = 0
    for g in genes:
        module = membership[g]
        if module != 'UNCLASSIFIED':
            if memberCount[module] < CML_ARGS.wgcnaParameters['minModuleSize']:
                membership[g] = 'UNCLASSIFIED'
            else:
                classifiedGeneCount = classifiedGeneCount + 1
                modules[module] = 1
    return membership, len(modules), classifiedGeneCount

def evaluate_fit(kME, membership, genes):
    '''
    evaluate eigengene similarity (KME) as
    a measure of goodness of fit
    label poorly fitting genes as unclassified
    then remove any modules whose n < minMoudleSize
    returns updated membership and number of modules
    and classified gene count (to test convergence
    conditions)
    '''
    memberCount = {}
    for g in genes:
        module = membership[g]

        if membership == 'UNCLASSIFIED':
            continue

        if module in memberCount:
            memberCount[module] = memberCount[module] + 1
        else:
            memberCount[module] = 1

        if kME[g] < CML_ARGS.wgcnaParameters['minKMEtoStay']:
            membership[g] = 'UNCLASSIFIED'
            memberCount[module] = memberCount[module] - 1

    membership, moduleCount, classifiedGeneCount = remove_small_modules(memberCount, membership,
                                                                        genes)
    return membership, moduleCount, classifiedGeneCount

def run_iteration(iteration, data, membership, kME):
    '''
    run iteration of blockwise WGCNA
    return membership and eigengenes
    '''

    # convergence criteria
    moduleCount = 0 # number of modules detected
    classifiedGeneCount = 0 # number of genes classified

    # run blockwise WGCNA
    blocks = run_wgcna(data, iteration)
    if CML_ARGS.saveBlocks:
        utils.saveObject(blocks, 'blocks', iteration + '-blocks.RData')

    # extract eigengenes from blocks
    # if there are eigengenes, then evaluate fitness
    eigengenes = utils.eigengenes(iteration, blocks, data.colnames)
    if eigengenes.nrow != 0:
        write_eigengenes(eigengenes)

        # extract module membership from blocks and update
        membership = update_membership(iteration, data.rownames, blocks, membership)
        write_membership(iteration, membership, False)

        # calculate eigengene connectivity (kME) to
        # assigned module
        kME = update_kme(kME, data, membership, eigengenes)
        write_kme(iteration, kME)

        # evaluate fit & update membership again
        # output pruned membership
        membership, moduleCount, classifiedGeneCount = evaluate_fit(kME, membership, data.rownames)
        write_membership(iteration, membership, True)

    return membership, kME, moduleCount, classifiedGeneCount

# ========================
# iWGCNA Main
# ========================

def iWGCNA():
    '''
    run iterative WGCNA
    '''
    # initialize iteration counters and iteration labels
    runId = 0
    passId = 0
    iteration = set_iteration_label(runId, passId)

    # initialize variables and data (to input DATA)
    kME = None
    membership = None
    passData = DATA # input data for pass
    iterationData = passData # input data for iteration

    # initialize convergence flags
    algConverged = False
    passConverged = False

    while not algConverged:
        iteration = set_iteration_label(runId, passId)

        # run an iteration of WGCNA + goodness of fit test
        membership, kME, moduleCount, classifiedGeneCount = run_iteration(
            iteration, iterationData, membership, kME)

        # if there are no residuals
        # (classified gene count = number of genes input)
        # then the pass has converged
        if classifiedGeneCount == iterationData.nrow:
            passConverged = True

        # if no modules were detected,
        # then the algorithm has converged
        if moduleCount == 0:
            algConverged = True
            logger.info("CLASSIFICATION COMPLETE")

        write_gene_counts(iteration, iterationData.nrow, classifiedGeneCount)

        if passConverged:
            log_pass_completion(passId, iteration, passData.nrow,
                                classifiedGeneCount, moduleCount)

            # set residuals of the pass as new
            # input dataset for the next pass
            passData = get_residuals(passData, membership)

            # and as the dataset for the current iteration
            iterationData = passData

            # increment pass id and reset run id
            passId = passId + 1
            runId = 0

            # reset pass convergence flag
            passConverged = False

        else:
            # other wise remove residuals and run again
            # with classified gene set to evaluate whether
            # classification can be improved after
            # connections to pruned genes are removed
            iterationData = remove_residuals(iterationData, membership)
            runId = runId + 1

    return membership, kME

def reassign_membership(membership, kME):
    '''
    Evaluate eigengene connectivity (kME)
    for each gene against each module
    eigengene.  If kME(module) > kME(assigned_module)
    and the p-value <= the reassignThreshold (of WGCNA
    parameters) then reassign the module
    membership of the gene
    '''

    eigengenes = read_data('eigengenes.txt')
    reassignmentCount = 0
    for module in eigengenes.rownames:
        # calculate kME of all genes to the module eigengene
        moduleEigengene = eigengenes.rx(module, True)
        moduleKME = calculate_kme(DATA, moduleEigengene, True)

        for gene in membership:
            currentModule = membership[gene]
            currentKME = kME[gene]

            newKME = moduleKME.rx2('cor').rx(gene, 1)[0]
            pvalue = moduleKME.rx2('p').rx(gene, 1)[0]

            # logger.debug('\t'.join((gene, currentModule, str(currentKME), module, str(newKME), str(pvalue))))

            if (currentModule == "UNCLASSIFIED" \
                and newKME >= CML_ARGS.wgcnaParameters['minKMEtoStay']) \
                or (newKME > currentKME \
                and pvalue < CML_ARGS.wgcnaParameters['reassignThreshold']):
                membership[gene] = module
                kME[gene] = newKME
                reassignmentCount = reassignmentCount + 1

    logger.info("Reassigned " + str(reassignmentCount) + " genes in final kME review.")
    return membership, kME

def main():
    '''
    main function:
    1. run iterative WGCNA
    2. make last check goodness of fit test to
       reassign module membership(s) after all
       iterations
    3. transpose output files so samples are in
       columns and genes in rows to improve
       readability
    '''

    # run iterative WGCNA
    membership, kME = iWGCNA()

    # use kME goodness of fit to reassign module
    # membership now that all modules are estimated
    membership, kME = reassign_membership(membership, kME)
    write_membership('final', membership, True)
    write_kme('final', kME)

    # transpose membership and kME files (so samples are columns)

if __name__ == '__main__':
    logger = None
    
    try:
        CML_ARGS = parse_command_line_args()

        # set up working environment
        # ===============================
        # import R packages
        wgcna = importr('WGCNA')
        igraph = importr('igraph')
        base = importr('base')
        stats = importr('stats')
        utils = SignatureTranslatedAnonymousPackage(rs.util_functions, 'utils')

        # create working directory and
        # initialize R workspace
        create_dir(CML_ARGS.workingDir)
        initialize_r_workspace()
    
        # initialize log
        logger = initialize_log()
        log_parameters()
        
        # load input data
        DATA = read_data(CML_ARGS.inputFile)
        log_input_data()

        main()

    finally:
        if logger is not None:
            logger.info('iWGCNA: FINISHED')
            logger.info(strftime("%c"))

__author__ = 'Emily Greenfest-Allen'
__copyright__ = 'Copyright 2016, University of Pennsylvania'
