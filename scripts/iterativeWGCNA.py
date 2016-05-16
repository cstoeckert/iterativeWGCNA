#!/usr/bin/env python2.7
"""
Perform iterative WGCNA analysis

python dependencies:
  * rpy2

R dependencies:
  * igraph
  * WGCNA
  * lattice
"""

import sys
from time import strftime
from iwgcna import membership, eigengenes, kme, expression
from iwgcna import warning, initialize_r_workspace, parse_command_line_args
from iwgcna import process as iwgcna
from iwgcna.utils import io as iwgcna_io
from iwgcna.utils import log

def iterative_wgcna():
    '''
    iteratively apply WGCNA until convergence
    conditions are met
    '''
  
    runId = 1
    passId = 1
  
    kmeMap = kme.initialize(DATA)
    membershipMap = membership.initialize(DATA)
    passData = DATA # input data for pass
    iterationData = passData # input data for iteration

    algConverged = False
    passConverged = False

    while not algConverged:
        iteration = iwgcna.set_iteration_label(runId, passId)

        if CML_ARGS.verbose:
            warning("Running iteration: " + iteration)

        # run an iteration of WGCNA + goodness of fit test
        membershipMap, kmeMap = iwgcna.run_iteration(
            iteration, iterationData, membershipMap, kmeMap,
            CML_ARGS.wgcnaParameters, CML_ARGS.saveBlocks)

        moduleCount = membership.count_modules(membershipMap)
        classifiedGeneCount = membership.count_classified_genes(membershipMap)

        # if there are no residuals
        # (classified gene count = number of genes input)
        # then the pass has converged
        if classifiedGeneCount == iterationData.nrow:
            passConverged = True

        # if no modules were detected,
        # then the algorithm has converged
        if moduleCount == 0:
            algConverged = True
            message = "No modules detected for iteration " + iteration \
                      + ". Classification complete."
            logger.info(message)
            if CML_ARGS.verbose:
                warning(message)

        iwgcna.write_gene_counts(iteration, iterationData.nrow, classifiedGeneCount)

        if passConverged:
            summary = log.pass_completion(logger, passId, iteration, passData.nrow,
                                          classifiedGeneCount, moduleCount)
            if CML_ARGS.verbose:
                warning(summary)

            # set residuals of the pass as new
            # input dataset for the next pass
            passData = expression.get_residuals(passData, membershipMap)

            # and as the dataset for the current iteration
            iterationData = passData

            # increment pass id and reset run id
            passId = passId + 1
            runId = 1

            # reset pass convergence flag
            passConverged = False

        else:
            # otherwise remove residuals and run again
            # with classified gene set to evaluate whether
            # classification can be improved after
            # connections to pruned genes are removed
            iterationData = expression.remove_residuals(iterationData, membershipMap)
            runId = runId + 1

    return membershipMap, kmeMap


def main():
    '''
    main function:
    1. run iterative WGCNA
    2. make last check goodness of fit test to
       reassign module membership(s) after all
       iterations
    3. transpose output files so samples are in
       columns and genes in rows to improve
       readability and so that files are EXCEL friendly
    '''

    if CML_ARGS.verbose:
        warning("Starting iWGCNA")
    # run iterative WGCNA
    membershipMap, kmeMap = iterative_wgcna()

    # use kME goodness of fit to reassign module
    # membership now that all modules are estimated
    if CML_ARGS.verbose:
        warning("Making final goodness of fit assessment")
    modules = membership.get_modules(membershipMap)
    logger.info("Found " + str(len(modules)) + " modules.")
    logger.info(modules)

    # extract final eigengene matrix
    eigengeneMatrix = eigengenes.load_from_file('eigengenes.txt')
    eigengeneMatrix = eigengenes.extract_modules(eigengeneMatrix, modules)
    logger.debug(eigengeneMatrix)
    eigengenes.write(eigengeneMatrix, True)

    reassignedCount, membershipMap, kmeMap = membership.best_fit(membershipMap, eigengeneMatrix, DATA,
                                                                 kmeMap, CML_ARGS.wgcnaParameters)
    logger.info("Reassigned " + str(reassignedCount) + " genes in final kME review.")
    if CML_ARGS.verbose:
        warning("Reassigned " + str(reassignedCount) + " genes in final kME review.")

    membership.write('final', membershipMap, True)
    kme.write('final', kmeMap)

    if CML_ARGS.verbose:
        warning("Generating final output")

    # transpose membership and kME files (so samples are columns)
    iwgcna_io.transpose_file_contents("pre-pruning-membership.txt")
    iwgcna_io.transpose_file_contents("membership.txt")
    iwgcna_io.transpose_file_contents("eigengene-connectivity.txt")

    logger.info('iWGCNA: SUCCESS')


if __name__ == '__main__':

    logger = None
    try:
        CML_ARGS = parse_command_line_args()

        if CML_ARGS.verbose:
            warning("Initializing workspace")

        # initialize R workspace and logs
        iwgcna_io.create_dir(CML_ARGS.workingDir)
        initialize_r_workspace(CML_ARGS.workingDir, CML_ARGS.allowWGCNAThreads)

        logger = log.initialize(CML_ARGS.workingDir)
        log.parameters(logger, CML_ARGS)

        # load input data
        if CML_ARGS.verbose:
            warning("Loading Data")

        # gives a weird R error that I'm having trouble catching
        # hence the nested block
        # TODO: identify the exact exception and remove
        # nested try-catch
        try:
            DATA = iwgcna_io.read_data(CML_ARGS.inputFile)
        except:
            logger.error("Unable to open input file: " + CML_ARGS.inputFile)
            sys.exit(1)

        log.input_data(logger, CML_ARGS.inputFile, DATA.ncol, DATA.nrow)

        main()

    except Exception:
        logger.exception('iWGCNA: FAILED')

    finally:
        if logger is not None:
            logger.info(strftime("%c"))
     

__author__ = 'Emily Greenfest-Allen'
__copyright__ = 'Copyright 2016, University of Pennsylvania'
