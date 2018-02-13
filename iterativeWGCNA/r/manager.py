# pylint: disable=invalid-name
# pylint: disable=no-self-use
'''
Wrappers for R functions;
performs conversions to rpy2.robjects
where required
'''

import logging
import rpy2.robjects as ro
from collections import OrderedDict

# import rpy2.rlike.container as rlc
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


    def heatmap_annotation_data_frame(self, categories, annotation):
        '''
        takes a dict of gene->value and creates a data frame
        data frame
        assume annotation is an ordered dict
        updates column names to names
        '''
        df = base().as_data_frame(base().t(ro.DataFrame(annotation)))
        df.colnames = ro.StrVector(categories)
      
        return df


    def heatmap_annotation_key(self, name, colors):
        '''
        generates data frame for color key for the annotation
        from a dict
        '''
        keyColors = ro.StrVector([c for c in colors.values()])
        keyColors.names = colors.keys()
        key = OrderedDict()
        key[name] = keyColors

        return ro.ListVector(key)


    def heatmap(self, clusterCols=False, params=None):
        '''
        plot a heatmap with options specified in params
        (see pheatmap documentation for all options)
        '''
        self.params['mat'] = base().as_matrix(self.log2())
        self.params['border'] = ro.NA_Logical
        self.params['cluster_cols'] = clusterCols
        if clusterCols:
            self.params['clustering_distance_cols'] = 'correlation'
        self.params['clustering_distance_rows'] = 'correlation'
        self.params['show_rownames'] = True if self.data.nrow <= 50 else False
        self.params['scale'] = 'row'
        self.params['color'] = rsnippets.BlWhRed()

        if params is not None:
            self.params.update(params)

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
        self.params['x'] = ro.FloatVector(self.data)
        self.params['labels'] = False

        if params is not None:
            self.params.update(params)

        graphics().hist(**self.params)

        if vline is not None:
            lineParams = {'v': vline, 'col': 'red'}
            graphics().abline(**lineParams)



