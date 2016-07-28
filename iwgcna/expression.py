'''
functions for manipulating expression matrices
'''

import rpy2.robjects as ro

class Expresson(object):
    def __init__(self, data):
        self.data = data
        self.size = len(self.data)
        return None

    def get_feature_expression(self, features):
        '''
        subsets expression data
        returning expression for list of features
        '''
        return self.data.rx(ro.StrVector(features), True)

    def get_residual_expression(self, unclassifiedFeatures):
        '''
        subsets expression data
        returning expression for only residuals to the fit
        '''
        if unclassifiedFeatures is None:
            return None
        else:
            return self.get_feature_expression(unclassifiedFeatures)

    def get_fit_expression(self, fitFeatures):
        '''
        subsets expression data
        return only features that passed the
        goodness of fit assessment
        '''
        if fitFeatures is None:
            return None
        else:
            return self.get_feature_expression(fitFeatures)


    def read_data(fileName):
        '''
        reads the data from a file
        '''
