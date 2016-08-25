# pylint: disable=invalid-name
# pylint: disable=no-self-use
'''
wgcna functions
'''

# import rpy2.robjects as ro
from .r.imports import base, wgcna

class WgcnaManager(object):
    '''
    wrappers for running WGCNA functions
    '''
    def __init__(self, data, params):
        self.exprData = data
        self.params = params
        self.adjacencyMatrix = None
        self.TOM = None
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


    def __transpose_data(self):
        '''
        transpose the data frame (required for some WGCNA functions)
        '''
        return base().t(self.exprData)


    def blockwise_modules(self):
        '''
        run blockwise WGCNA
        '''
        self.params['datExpr'] = self.__transpose_data()
        blocks = wgcna().blockwiseModules(**self.params)
        self.collect_garbage()
        return blocks


    def collect_garbage(self):
        '''
        run WGCNA garbage collection
        '''
        wgcna().collectGarbage()


    def adjacency(self):
        '''
        calculate adjacency matrix; from pearson correlation
        '''
        adjParams = {}
        adjParams['power'] = self.params['power'] if 'power' in self.params else 6
        adjParams['corFunc'] = 'cor'
        adjParams['corOptions'] = "use='p'"
        adjParams['exprData'] = self.__transpose_data()

        self.adjacencyMatrix = wgcna().adjacency(**adjParams)
        self.collect_garbage()


    def TOM_dist(self):
        '''
        calculate dis-Topological Overlap Matrix
        '''
        self.TOM = wgcna().TOMdist(self.adjacencyMatrix)
        self.collect_garbage()

