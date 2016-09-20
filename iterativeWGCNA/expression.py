# pylint: disable=invalid-name

'''
functions for manipulating expression profile matrices
'''

import rpy2.robjects as ro

class Expression(object):
    '''
    store and manipulate expression profile matrix
    '''
    def __init__(self, data):
        self.profiles = data
        self.size = len(self.profiles)
        return None


    def genes(self):
        '''
        return genes (row names)
        '''
        return self.profiles.rownames


    def nrow(self):
        '''
        return number of rows
        '''
        return self.profiles.nrow


    def ncol(self):
        '''
        return number of columns
        '''
        return self.profiles.ncol


    def samples(self):
        '''
        return column names (samples)
        '''
        return self.profiles.colnames


    def expression(self):
        '''
        wrapper for accessing self.profiles
        '''
        return self.profiles


    def gene_expression(self, genes):
        '''
        subsets expression data
        returning expression for list of genes
        '''
        return self.profiles.rx(ro.StrVector(genes), True)


    def residual_expression(self, unclassifiedGenes):
        '''
        subsets expression data
        returning expression for only residuals to the fit
        '''
        if unclassifiedGenes is None:
            return None
        else:
            return self.gene_expression(unclassifiedGenes)


    def fit_expression(self, fitGenes):
        '''
        subsets expression data
        return only genes that passed the
        goodness of fit assessment
        '''
        if fitGenes is None:
            return None
        else:
            return self.gene_expression(fitGenes)
