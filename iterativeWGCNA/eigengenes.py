# pylint: disable=invalid-name
# pylint: disable=unused-import
'''
manage eigengenes
'''

from __future__ import print_function

import logging

import rpy2.robjects as ro
from .r.imports import base, stats, rsnippets
from .io.utils import write_data_frame
from .wgcna import WgcnaManager

class Eigengenes(object):
    '''
    manage and manipulate eigengene matrices
    '''

    def __init__(self, matrix=None, debug=False):
        self.debug = debug
        self.logger = logging.getLogger('iterativeWGCNA.Eigengenes')
        self.matrix = matrix


    def extract_from_blocks(self, iteration, blocks, samples):
        '''
        extract eigenenges from blockwise WGCNA results
        '''
        self.matrix = rsnippets.extractEigengenes(iteration, blocks, samples)


    def samples(self):
        '''
        return sample names
        '''
        return self.matrix.names


    def nrows(self):
        '''
        wrapper for returning number of rows in
        the eigengene matrix
        '''

        return self.matrix.nrow


    def load_matrix_from_file(self, fileName):
        '''
        loads eigengenes from file into an R DataFrame
        '''
        self.matrix = ro.DataFrame.from_csvfile(fileName, sep='\t',
                                                header=True, row_names=1)


    def write(self, prefix=''):
        '''
        writes the eigengene matrix to file
        '''
        fileName = prefix + 'eigengenes.txt'
        write_data_frame(self.matrix, fileName, 'Module')


    def similarity(self, module=None):
        '''
        calculate similarity between eigengene for a specific
        module and all the other eigengenes

        if no module is specified, calculate the similarity matrix
        between all eigengenes
        '''
        if module is None:
            sim = base().as_data_frame(stats().cor(base().t(self.matrix)))
        else:
            sim = base().as_data_frame(stats().cor(base().t(self.matrix), \
                                            base().t(self.matrix.rx(module, True))))
        return sim


    def correlation(self, m1, m2):
        '''
        calculate correlation between two module eigengenes
        '''
        e1 = self.get_module_eigengene(m1)
        e2 = self.get_module_eigengene(m2)
        cor = base().as_data_frame(stats().cor(base().t(e1), base().t(e2)))
        cor = round(cor.rx(1, 1)[0], 1)
        return cor


    def equal(self, m1, m2, threshold=0.0):
        '''
        check if 2 module eigengenes are "equivalent"
        (1 - correlation <= threshold)
        '''
        cor = self.correlation(m1, m2)
        return 1.0 - cor <= threshold


    def get_module_eigengene(self, module):
        '''
        return a module eigengene
        '''
        return self.matrix.rx(module, True)


    def extract_subset(self, modules):
        '''
        return a submatrix
        '''
        if self.debug:
            self.logger.debug("Extracting eigengenes for the following modules:")
            self.logger.debug(modules)

        if self.debug:
            self.logger.debug("Converting module list to ro.StrVector; see R-log")
            ro.r("print('Converting module list to ro.StrVector to extract eigengenes:')")

        vector = ro.StrVector(modules)

        if self.debug:
            self.logger.debug(vector)

        if self.debug:
            self.logger.debug("Extracted submatrix, see R-log")
            ro.r("print('Extracted eigengene submatrix:')")


        newMatrix = self.matrix.rx(vector, True)

        if self.debug:
            self.logger.debug(newMatrix)

        return newMatrix



    def is_empty(self):
        '''
        return True if matrix is empty
        '''
        return self.matrix.nrow == 0


    def update_to_subset(self, modules):
        '''
        update matrix to subset specified by modules
        '''
        self.matrix = self.extract_subset(modules)


    def recalculate(self, profiles, membership, power=6):
        '''
        recalculate eigengenes given membership
        and profiles
        '''
        manager = WgcnaManager(profiles, {'power':power}, debug=self.debug)

        self.matrix = rsnippets.extractRecalculatedEigengenes(
            manager.module_eigengenes(membership.values()),
            self.samples())
