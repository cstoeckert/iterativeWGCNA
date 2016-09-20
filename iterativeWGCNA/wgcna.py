# pylint: disable=invalid-name
# pylint: disable=no-self-use
'''
wgcna functions
'''

import logging
import rpy2.robjects as ro
from .r.imports import base, wgcna, rsnippets, stats

class WgcnaManager(object):
    '''
    wrappers for running WGCNA functions
    '''
    def __init__(self, data, params):
        self.data = data
        self.params = params
        self.adjacencyMatrix = None
        self.TOM = None
        self.dissimilarityMatrix = None
        self.geneTree = None
        self.logger = logging.getLogger('iterativeWGCNA.WgcnaManager')
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
        return base().t(self.data)


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


    def adjacency(self, signed=True, removeNegatives=False, removeSelfReferences=False):
        '''
        calculate adjacency matrix; from pearson correlation
        '''
        adjParams = {}
        adjParams['power'] = self.params['power'] if 'power' in self.params else 6
        adjParams['corFnc'] = 'cor'
        adjParams['corOptions'] = "use='p'"
        if signed:
            adjParams['type'] = 'signed'
        else:
            adjParams['type'] = 'unsigned'
        adjParams['datExpr'] = self.__transpose_data()

        self.adjacencyMatrix = wgcna().adjacency(**adjParams)
        self.collect_garbage()
        if removeNegatives:
            self.adjacencyMatrix = rsnippets.filterNegatives(self.adjacencyMatrix)

        if removeSelfReferences:
            self.adjacencyMatrix = rsnippets.diag(self.adjacencyMatrix, 0)


    def TOM_dist(self):
        '''
        calculate Topological Overlap Matrix from adjacency matrix
        '''
        self.TOM = wgcna().TOMdist(self.adjacencyMatrix)
        self.collect_garbage()


    def TOM_similarity_from_expr(self):
        '''
        calculate Topological Overlap Matrix from expression data
        '''
        funcParams = {}
        funcParams['power'] = self.params['power'] if 'power' in self.params else 6
        funcParams['datExpr'] = self.__transpose_data()
        self.TOM = wgcna().TOMsimilarityFromExpr(**funcParams)


    def plot_network_heatmap(self, genes, title="Network Heatmap", useTOM=False):
        '''
        wrapper for plotNetworkHeatmap function
        by default plots correlation, set useTOM to true
        to plot topological overlap instead
        '''

        funcParams = {}
        funcParams['power'] = self.params['power'] if 'power' in self.params else 6
        funcParams['networkType'] = self.params['networkType'] \
          if 'networkType' in self.params else 'unsigned'
        funcParams['main'] = title
        funcParams['datExpr'] = self.__transpose_data()
        funcParams['plotGenes'] = ro.StrVector(genes)
        funcParams['useTOM'] = useTOM

        wgcna().plotNetworkHeatmap(**funcParams)


    def generate_gene_tree(self):
        '''
        generate hierarchical cluster from dissimilarity matrix
        '''
        distMatrix = stats().as_dist(self.dissimilarityMatrix)
        self.geneTree = stats().hclust(distMatrix, method="average")


    def plot_network_overview(self, moduleColors, title, useTOM=False):
        '''
        wrapper for TOMplot which provides a graphical representation
        of the Topological Overlap Matrix or correlation matrix using a heatmap
        and hierarchical clustering dendrogram annotated by
        module colors
        '''

        if useTOM:
            self.TOM_similarity_from_expr()
            self.dissimilarityMatrix = rsnippets.dissMatrix(self.TOM)
            self.dissimilarityMatrix = rsnippets.powerWeightMatrix(self.TOM, 7)
        else:
            self.adjacency()
            self.dissimilarityMatrix = rsnippets.dissMatrix(self.adjacencyMatrix)

        # self.dissimilarityMatrix = rsnippets.diag(self.dissimilarityMatrix, ro.NA_Integer)
        self.generate_gene_tree()

        funcParams = {}
        funcParams['dissim'] = self.dissimilarityMatrix
        funcParams['dendro'] = self.geneTree
        funcParams['Colors'] = moduleColors
        funcParams['main'] = title
        wgcna().TOMplot(**funcParams)


    def plot_eigengene_network(self):
        '''
        wrapper for plotEigengeneNetworks function
        plots an eigengene network
        '''
        funcParams = {}
        funcParams['multiME'] = base().as_data_frame(self.__transpose_data())
        funcParams['setLabels'] = ''
        funcParams['marDendro'] = ro.IntVector([0, 4, 1, 2])
        funcParams['marHeatmap'] = ro.IntVector([3, 4, 1, 2])
        funcParams['cex.lab'] = 0.8
        funcParams['xLabelsAngle'] = 90
        funcParams['colorLabels'] = False
        funcParams['signed'] = True

        wgcna().plotEigengeneNetworks(**funcParams)

