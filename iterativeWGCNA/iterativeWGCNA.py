# pylint: disable=invalid-name
# pylint: disable=bare-except
# pylint: disable=broad-except
# pylint: disable=too-many-instance-attributes
'''
main application
'''

from __future__ import print_function

import logging
import sys
import os
from time import strftime

import rpy2.robjects as ro
from .genes import Genes
from .expression import Expression
from .eigengenes import Eigengenes
from .network import Network
from .wgcna import WgcnaManager
from .io.utils import create_dir, read_data, warning, write_data_frame
from .r.imports import base, wgcna, rsnippets


class IterativeWGCNA(object):
    '''
    main application

    flag report = True when generating
    result from existing output
    '''

    def __init__(self, args, report=False):
        self.args = args
        create_dir(self.args.workingDir)
        if not report:
            self.__verify_clean_working_dir()

        if report == 'merge':
            self.args.enableWGCNAThreads = False

        self.__initialize_log(report)
        self.logger.info(strftime("%c"))

        self.__initialize_R(report)
        if not report:
            self.__log_parameters()

        if self.args.debug:
            warning("Running in DEBUG mode.")
            warning("Rpy2 will print debugging messages and variable (e.g., matrix/vector) contents to the R log")
            warning("Thus, empty debug statements 'DEBUG:      ' in the iterativeWGCNA log should have a corresponding output in the R log")

        # load expression data and
        # initialize Genes object
        # to store results
        self.profiles = None
        self.__load_expression_profiles()
        self.__log_input_data()
        self.genes = Genes(self.profiles, debug=self.args.debug)
        self.eigengenes = Eigengenes(debug=args.debug)
        self.modules = None # will be hash of module name to color for plotting

        if not report:
            self.passCount = 1
            self.iterationCount = 1
            self.iteration = None # unique label for iteration
            self.algorithmConverged = False
            self.passConverged = False


    def __verify_clean_working_dir(self):
        '''
        verifies that working directory does not contain
        iterativeWGCNA output files
        exits to avoid accidental overwrite of earlier runs
        '''
        conflictingFiles = set(('final-eigengenes.txt', 'final-membership.txt'))
        files = set(os.listdir(self.args.workingDir))
        if len(files.intersection(conflictingFiles)) > 0:
            warning("Working Directory: " + self.args.workingDir \
                               + " contains final output from a prior run of iterativeWGCNA.  Exiting...")
            sys.exit(1)


    def run_pass(self, passGenes):
        '''
        run a single pass of iterative WGCNA
        (prune data until no more residuals are found)
        '''

        passDirectory = 'pass' + str(self.passCount)
        create_dir(passDirectory)
        write_data_frame(self.profiles.gene_expression(passGenes),
                         os.path.join(passDirectory, 'initial-pass-expression-set.txt'),
                         'Gene')

        iterationGenes = passGenes

        while not self.passConverged:
            self.run_iteration(iterationGenes)

            moduleCount = self.genes.count_modules(iterationGenes)
            classifiedGeneCount = self.genes.count_classified_genes(iterationGenes)

            self.write_run_summary(len(iterationGenes), classifiedGeneCount)

            # if there are no residuals
            # (classified gene count = number of genes input)
            # then the pass has converged
            if classifiedGeneCount == len(iterationGenes):
                self.passConverged = True
                self.__summarize_classification(passDirectory + '/')
            else:
                # run again with genes classified in current pass
                iterationGenes = self.genes.get_classified_genes(iterationGenes)
                self.iterationCount = self.iterationCount + 1

            # if no modules were detected,
            # then the algorithm has converged
            # exit the pass
            if moduleCount == 0:
                self.algorithmConverged = True
                self.passConverged = True
                self.__log_alogorithm_converged()


    def run_iterative_wgcna(self):
        '''
        run iterative WGCNA
        '''
        if self.args.verbose:
            warning("Beginning iterations")

        # genes involved in current iteration
        passGenes = self.profiles.genes()

        while not self.algorithmConverged:
            self.run_pass(passGenes)
            classifiedGeneCount = self.genes.count_classified_genes(passGenes)
            self.__log_pass_completion()
            self.__log_gene_counts(len(passGenes), classifiedGeneCount)

            if not self.algorithmConverged:
                # set residuals as new gene list
                passGenes = self.genes.get_unclassified_genes()

                # increment pass counter and reset iteration counter
                self.passCount = self.passCount + 1
                self.iterationCount = 1

                # reset pass convergence flag
                self.passConverged = False

        self.iteration = 'FINAL'
        self.genes.iteration = self.iteration
        self.__log_gene_counts(self.genes.size, self.genes.count_classified_genes())
        self.__summarize_classification('final-')

        # output current eigengenes for all modules, not just ones from last pass
        self.eigengenes.load_matrix_from_file('eigengenes.txt')
        modules = self.genes.get_modules()
        self.eigengenes.update_to_subset(modules)
        self.eigengenes.write('final-')

        self.iteration = 'MERGED'
        self.genes.iteration = self.iteration
        self.merge_close_modules()
        self.reassign_genes_to_best_fit_module()

        self.__log_gene_counts(self.genes.size, self.genes.count_classified_genes())

        self.__summarize_classification('merged-' + str(self.args.finalMergeCutHeight) + '-')
        self.eigengenes.write('merged-' + str(self.args.finalMergeCutHeight) + '-')
        os.remove("eigengenes.txt")


    def merge_close_modules_from_output(self):
        '''
        load data from output and remerge
        '''
        self.genes.load_membership()
        self.merge_close_modules('final-')
        self.reassign_genes_to_best_fit_module()
        self.__log_gene_counts(self.genes.size, self.genes.count_classified_genes())
        self.genes.write('adjusted-merge-' + str(self.args.finalMergeCutHeight) + '-')
        self.eigengenes.write('adjusted-merge-' + str(self.args.finalMergeCutHeight) + '-')
        # self.transpose_output_files()


    def summarize_results(self):
        '''
        generate summary output and graphics
        '''
        network = Network(self.args)
        network.build(self.genes, self.eigengenes)
        network.summarize_network()


    def run(self):
        '''
        main function --> makes calls to run iterativeWGCNA,
        catches errors, and logs time
        '''

        try:
            self.run_iterative_wgcna()
            # self.summarize_results() # can cause memory issues so, removing
            self.logger.info('iterativeWGCNA: SUCCESS')
        except Exception:
            if self.logger is not None:
                self.logger.exception('iterativeWGCNA: FAIL')
            else:
                raise
        finally:
            if self.logger is not None:
                self.logger.info(strftime("%c"))


    def reassign_genes_to_best_fit_module(self):
        '''
        use kME goodness of fit to reassign module
        membership
        '''
        if self.args.verbose:
            warning("Making final goodness of fit assessment")

        count = self.genes.reassign_to_best_fit(self.eigengenes,
                                                self.args.wgcnaParameters['reassignThreshold'],
                                                self.args.wgcnaParameters['minKMEtoStay'])
        self.logger.info("Reassigned " + str(count) + " genes in final kME review.")
        if self.args.verbose:
            warning("Reassigned " + str(count) + " genes in final kME review.")


    def merge_close_modules(self, prefix=''):
        '''
        merge close modules based on similiarity in eigengenes
        update membership, kME, and eigengenes accordingly
        '''
        if self.args.verbose:
            warning("Extracting final eigengenes and merging close modules")

        modules = self.genes.get_modules()
        self.__log_final_modules(modules)

        self.eigengenes.load_matrix_from_file(prefix + 'eigengenes.txt')
        self.eigengenes.update_to_subset(modules)

        self.eigengenes = self.genes.merge_close_modules(self.eigengenes,
                                                         self.args.finalMergeCutHeight)


    def run_iteration(self, iterationGenes):
        '''
        run an iteration of blockwise WGCNA
        '''
        self.__generate_iteration_label()

        iterationDir = os.path.join('pass' + str(self.passCount), 'i' + str(self.iterationCount))
        create_dir(iterationDir)

        if self.args.verbose:
            warning("Iteration: " + self.iteration)

        self.genes.iteration = self.iteration
        iterationProfiles = self.profiles.gene_expression(iterationGenes)

        blocks = self.run_blockwise_wgcna(iterationProfiles, iterationDir)
        if not self.args.skipSaveBlocks:
            rsnippets.saveBlockResult(blocks, iterationProfiles,
                                      os.path.join(iterationDir, 'wgcna-blocks.RData'))

        # update eigengenes from blockwise result
        # if eigengenes are present (modules detected), evaluate
        # fitness and update gene module membership
        self.eigengenes.extract_from_blocks(self.iteration, blocks,
                                            self.profiles.samples())

        if not self.eigengenes.is_empty():
            self.eigengenes.write() # need to keep single file across all iterations
            self.eigengenes.write(iterationDir + '/')

            # extract membership from blocks and calc eigengene connectivity
            self.genes.update_membership(iterationGenes, blocks)
            self.genes.update_kME(self.eigengenes, iterationGenes)
            self.__summarize_classification(os.path.join(iterationDir, 'wgcna-'))

            self.genes.evaluate_fit(self.args.wgcnaParameters['minKMEtoStay'],
                                    iterationGenes)
            self.genes.remove_small_modules(self.args.wgcnaParameters['minModuleSize'])
            self.__summarize_classification(os.path.join(iterationDir, ''), True)


    def __summarize_classification(self, prefix, inclCounts=False):
        '''
        output gene summaries for the iteration
        incl: text summary of iteration, updated gene membership, kme histogram
        '''
        if 'final' in prefix or 'merge' in prefix:
            self.genes.write(prefix)
        else:
            self.genes.write(prefix, self.iteration)
        self.genes.plot_kme_histogram(self.iteration, prefix,
                                      self.args.wgcnaParameters['minKMEtoStay'])
        if inclCounts:
            self.genes.write_iteration_counts(prefix)


    def run_blockwise_wgcna(self, exprData, workingDir):
        '''
        run WGCNA
        '''
        manager = WgcnaManager(exprData, self.args.wgcnaParameters)
        manager.set_parameter('saveTOMFileBase', os.path.join(workingDir, self.iteration + '-TOM'))
        return manager.blockwise_modules()


    def __generate_iteration_label(self):
        '''
        generates the unique label for the iteration
        '''
        self.iteration = 'P' + str(self.passCount) + '_I' + str(self.iterationCount)


    def __load_expression_profiles(self):
        # gives a weird R error that I'm having trouble catching
        # when it fails
        # TODO: identify the exact exception
        try:
            self.profiles = Expression(read_data(self.args.inputFile))
        except:
            self.logger.error("Unable to open input file: " + self.args.inputFile)
            sys.exit(1)


    def __initialize_R(self, logType='run'):
        '''
        initialize R workspace and logs
        '''
        # set working directory
        base().setwd(self.args.workingDir)

        # suppress warnings
        ro.r['options'](warn=-1)

        # r log
        logFile = 'iterativeWGCNA-R.log'
        if logType == 'merge':
            logFile = 'adjust-merge-' + str(self.args.finalMergeCutHeight) + '-' + logFile

        rLogger = base().file(logFile, open='wt')
        base().sink(rLogger, type=base().c('output', 'message'))

        if self.args.enableWGCNAThreads:
            wgcna().enableWGCNAThreads()


    def __initialize_log(self, logType='run'):
        '''
        initialize log by setting path and file format
        '''
        logName = 'iterativeWGCNA.log'
        if logType == 'summary':
            logName = 'summarize-network-' + logName
        elif logType == 'merge':
            logName = 'adjust-merge-' + str(self.args.finalMergeCutHeight) + '-' + logName

        logging.basicConfig(filename=self.args.workingDir + '/' + logName,
                            filemode='w', format='%(levelname)s: %(message)s',
                            level=logging.DEBUG)

        logging.captureWarnings(True)
        self.logger = logging.getLogger(__name__)


    def __log_alogorithm_converged(self):
        '''
        log algorithm convergence
        '''
        message = "No modules detected for iteration " + self.iteration \
                  + ". Classification complete."
        self.logger.info(message)
        if self.args.verbose:
            warning(message)


    def __log_parameters(self):
        '''
        log WGCNA parameter choices and working
        directory name
        '''

        self.logger.info("Working directory: " + self.args.workingDir)
        self.logger.info("Saving blocks for each iteration? "
                         + ("FALSE" if self.args.skipSaveBlocks else "TRUE"))
        self.logger.info("Merging final modules if cutHeight <= " + str(self.args.finalMergeCutHeight))
        self.logger.info("Allowing WGCNA Threads? "
                         + ("TRUE" if self.args.enableWGCNAThreads else "FALSE"))
        self.logger.info("Running WGCNA with the following params:")
        self.logger.info(self.args.wgcnaParameters)

        if self.args.verbose:
            warning("Working directory: " + self.args.workingDir)
            warning("Allowing WGCNA Threads? "
                    + ("TRUE" if self.args.enableWGCNAThreads else "FALSE"))
            warning("Merging final modules if cutHeight <= " + str(self.args.finalMergeCutHeight))
            warning("Running WGCNA with the following params:")
            warning(self.args.wgcnaParameters)


    def __log_input_data(self):
        '''
        log input details
        '''
        self.logger.info("Loaded file: " + self.args.inputFile)
        self.logger.info(str(self.profiles.ncol()) + " Samples")
        self.logger.info(str(self.profiles.nrow()) + " Genes")
        if self.args.verbose:
            warning("Loaded file: " + self.args.inputFile)
            warning(str(self.profiles.ncol()) + " Samples")
            warning(str(self.profiles.nrow()) + " Genes")


    def __log_pass_completion(self):
        '''
        summarize pass in log when convergence condition is met
        '''
        message = "Pass " + str(self.passCount) + " converged on iteration " \
                  + str(self.iterationCount) + "."
        self.logger.info(message)
        if self.args.verbose:
            warning(message)


    def __log_gene_counts(self, initialGeneCount, classifiedGeneCount):
        '''
        log classified gene count
        '''
        message = "FIT: " + str(classifiedGeneCount) + "; RESIDUAL: "
        message = message + str(initialGeneCount - classifiedGeneCount)
        self.logger.info(message)
        if self.args.verbose:
            warning(message)


    def __log_final_modules(self, modules):
        '''
        log modules
        '''
        message = "Found " + str(len(modules)) + " modules."
        self.logger.info(message)
        self.logger.info(modules)
        if self.args.verbose:
            warning(message)
            warning(modules)


    def write_run_summary(self, initial, fit):
        '''
        writes the number of kept and dropped genes at the end of an iteration
        '''
        fileName = 'iterative-wgcna-run-summary.txt'
        try:
            os.stat(fileName)
        except OSError:
            header = ('Iteration', 'Initial', 'Fit', 'Residual')
            with open(fileName, 'a') as f:
                print('\t'.join(header), file=f)
        finally:
            with open(fileName, 'a') as f:
                print('\t'.join((self.iteration, str(initial),
                                 str(fit), str(initial - fit))), file=f)
