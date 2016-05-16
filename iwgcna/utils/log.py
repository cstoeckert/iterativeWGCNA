'''
functions for initialing and handling logging
'''

import logging
from time import strftime
import os
from .io import warning

def initialize(workingDirectory):
    '''
    initialize log by setting path and file format
    '''
    logging.basicConfig(filename=os.path.join(workingDirectory, 'iwgcna.log'),
                        filemode='w', format='%(levelname)s: %(message)s',
                        level=logging.DEBUG)

    logging.captureWarnings(True)

    return logging.getLogger(__name__)


def parameters(params):
    '''
    log WGCNA parameter choices and working directory name
    '''
    logging.info(strftime("%c"))
    logging.info("Working directory: " + params.workingDir)
    logging.info("Allowing WGCNA Threads? "
                + ("TRUE" if params.allowWGCNAThreads else "FALSE"))
    logging.info("Running WGCNA with the following params:")
    logging.info(params.wgcnaParameters)

    if params.verbose:
        warning("Working directory: " + params.workingDir)
        warning("Allowing WGCNA Threads? "
                + ("TRUE" if params.allowWGCNAThreads else "FALSE"))
        warning("Running WGCNA with the following params:")
        warning(params.wgcnaParameters)


def input_data(fileName, nSamples, nGenes):
    '''
    log size of input data
    '''
    logging.info("Loaded file: " + fileName)
    logging.info(str(nSamples) + " Samples")
    logging.info(str(nGenes) + " Genes")


def pass_completion(passId, iteration, initialGeneCount,
                    classifiedGeneCount, moduleCount):
    '''
    summarize pass in log when convergence condition is met
    '''
    message = "Pass " + str(passId) + " converged on iteration " + iteration + "."
    message = message + " FIT: " + str(classifiedGeneCount) + "; RESIDUAL: "
    message = message + str(initialGeneCount - classifiedGeneCount)
    message = message + "; MODULES: " + str(moduleCount)
    logging.info(message)

    return message
