# pylint: disable=invalid-name
# pylint: disable=no-self-use
'''
R functions
'''

import logging
import rpy2.robjects as ro
from .imports import base, pheatmap, graphics, rsnippets

class RManager(object):
    '''
    wrappers for running R functions
    '''
    def __init__(self, data, params=None):
        self.logger = logging.getLogger('iterativeWGCNA.RManager')
        self.data = data
        if params is None:
            self.params = {}
        else:
            self.params = params

        return None


    def update_parameters(self, params):
        '''
        update/replace all parameters
        '''
        self.params = params


    def set_parameter(self, name, value):
        '''
        add or update a single parameter
        '''
        self.params[name] = value


    def remove_parameter(self, name):
        '''
        remove named parameter
        '''
        del self.params[name]


    def transpose_data(self):
        '''
        transpose the data frame (required for some WGCNA functions)
        '''
        return base().t(self.data)


    def log2(self):
        '''
        log2 data
        '''
        return base().log2(rsnippets.add(self.data, 1))


    def row_names(self):
        '''
        wrapper for getting
        row names (usually genes)
        '''
        return self.data.rownames


    def col_names(self):
        '''
        wrapper for getting column names
        (usually samples)
        '''
        return self.data.names


    def plot_heatmap(self, clusterCols=False, params=None):
        '''
        plot a heatmap with options specified in params
        (see pheatmap documentation for all options)
        '''
        if params is not None:
            self.params.update(params)

        self.params['mat'] = base().as_matrix(self.log2())
        self.params['border'] = ro.NA_Logical
        self.params['cluster_cols'] = clusterCols
        if clusterCols:
            self.params['clustering_distance_cols'] = 'correlation'
        self.params['clustering_distance_rows'] = 'correlation'
        self.params['show_rownames'] = True if self.data.nrow <= 50 else False
        self.params['scale'] = 'row'
        pheatmap().pheatmap(**self.params)


    def barchart(self, params=None):
        '''
        barchart
        '''

        self.params['height'] = base().as_numeric(self.data)
        if params is not None:
            self.params.update(params)

        graphics().barplot(**self.params)


    def histogram(self, vline=None, params=None):
        '''
        plot histogram with vline at x=vline
        '''
        if params is not None:
            self.params.update(params)

        self.logger.debug(self.data)
        graphics().hist(ro.FloatVector(self.data), main='test')
#        graphics().hist(ro.FloatVector(self.data), breaks=base().seq(0,1,0.1),
 #                       **self.params)

        if vline is not None:
            lineParams = {'v': vline, 'col': 'red'}
            graphics().abline(**lineParams)



