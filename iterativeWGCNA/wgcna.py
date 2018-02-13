# pylint: disable=invalid-name
# pylint: disable=no-self-use
'''
wgcna functions
'''

import logging
from collections import OrderedDict
import rpy2.robjects as ro
from .r.imports import base, wgcna, rsnippets, stats
from .r.manager import RManager

class WgcnaManager(RManager):
    '''
    wrappers for running WGCNA functions
    an extension of the RManager
    '''
    def __init__(self, data, params, debug=False):
        RManager.__init__(self, data, params)
        self.logger = logging.getLogger('iterativeWGCNA.WgcnaManager')
        self.adjacencyMatrix = None
        self.TOM = None
        self.dissimilarityMatrix = None
        self.geneTree = None
        self.moduleColors = None
        self.debug = debug
        return None


    def set_module_colors(self, modules):
        '''
        set module colors from dict of module properties
        expects dict to have a 'color' sub-dict
        '''
        self.moduleColors = OrderedDict((module, values['color']) \
                                            for module, values in modules.items())


    def blockwise_modules(self):
        '''
        run blockwise WGCNA
        '''
        self.params['datExpr'] = self.transpose_data()
        blocks = wgcna().blockwiseModules(**self.params)
        self.collect_garbage()
        return blocks


    def collect_garbage(self):
        '''
        run WGCNA garbage collection
        '''
        wgcna().collectGarbage()


    def adjacency(self, networkType='signed', removeNegatives=False, removeSelfReferences=False):
        '''
        calculate adjacency matrix; from pearson correlation
        '''

        adjParams = {}
        adjParams['power'] = self.params['power'] if 'power' in self.params else 6
        adjParams['corFnc'] = 'cor'
        adjParams['corOptions'] = "use='p'"
        adjParams['type'] = networkType
        adjParams['datExpr'] = self.transpose_data()

        self.adjacencyMatrix = wgcna().adjacency(**adjParams)
        self.collect_garbage()
        if removeNegatives:
            self.adjacencyMatrix = rsnippets.filterByThreshold(self.adjacencyMatrix, 0)

        if removeSelfReferences:
            self.adjacencyMatrix = rsnippets.diag(self.adjacencyMatrix, 0)

        self.collect_garbage()


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
        funcParams['datExpr'] = self.transpose_data()
        self.TOM = wgcna().TOMsimilarityFromExpr(**funcParams)
        self.collect_garbage()


    def plot_network_heatmap(self, membership, title, useTOM=False):
        '''
        plot network heatmap
        recapitulates WGCNA plotNetworkOverview to work around
        cutree issues
        uses pheatmap
        '''

        # TODO fix useTOM option for drawing with pheatmap and use of similarity
        # instead of dissimilarity

        # if useTOM:
            # self.TOM_similarity_from_expr()
            # self.dissimilarityMatrix = rsnippets.dissMatrix(self.TOM)
			# raising TOM to power of 7 recommended by WGCNA documentation
            # self.dissimilarityMatrix = rsnippets.powerWeightMatrix(self.TOM, 7)
        #else:

        self.adjacency()

        annotation = self.heatmap_annotation_data_frame(['Module'], membership)
        annotationKey = self.heatmap_annotation_key('Module', self.moduleColors)
        # manager = RManager(self.dissimilarityMatrix, None)
        self.heatmap(clusterCols=True, params={'scale': 'none',
                                               'mat': self.adjacencyMatrix,
                                               'show_colnames': False,
                                               'color': rsnippets.WhYlRed(),
                                               'annotation_col': annotation,
                                               'annotation_row': annotation,
                                               'annotation_colors': annotationKey})


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

        params = {}
        params['dissim'] = self.dissimilarityMatrix
        params['dendro'] = self.geneTree
        params['Colors'] = moduleColors
        params['main'] = title
        wgcna().TOMplot(**params)


    def module_eigengenes(self, membership):
        '''
        wrapper for moduleEigengenes function
        calculates eigengenes from profiles &
        module membership (gene -> membership dict)
        '''

        if self.debug:
            self.logger.debug("Running WGCNA moduleEigengenes function")
            self.logger.debug("Module assignments:")
            self.logger.debug(membership)

        params = {}
        params['softPower'] = self.params['power'] if 'power' in self.params else 6
        params['expr'] = base().as_data_frame(self.transpose_data())

        if self.debug:
            self.logger.debug("Converting membership list to ro.StrVector; see R-log")
            ro.r("print('Converting membership list to ro.StrVector for WGCNA moduleEigengenes:')")

        params['colors'] = ro.StrVector(list(membership))

        if self.debug:
            self.logger.debug(params['colors'])

        return wgcna().moduleEigengenes(**params)


    def plot_eigengene_network(self):
        '''
        wrapper for plotEigengeneNetworks function
        plots an eigengene network
        '''
        params = {}
        params['multiME'] = base().as_data_frame(self.transpose_data())
        params['setLabels'] = ''
        params['marDendro'] = ro.IntVector([0, 4, 1, 2])
        params['marHeatmap'] = ro.IntVector([3, 4, 1, 2])
        params['cex.lab'] = 0.8
        params['xLabelsAngle'] = 90
        params['colorLabels'] = False
        params['signed'] = True

        wgcna().plotEigengeneNetworks(**params)

