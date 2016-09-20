# pylint: disable=invalid-name

'''
manage network (for summaries)
'''

import logging
import sys
# import igraph
# import pandas as pd
# from rpy2.robjects import pandas2ri
from random import randint

from .wgcna import WgcnaManager
from .r.imports import grdevices

# from .eigengenes import Eigengenes
# from .genes import Genes
# from .module import Module

class Network(object):
    '''
    a network contains the classified genes
    and their assignments
    hash of modules and properties
    for summary purposes
    '''

    def __init__(self, genes, eigengenes, args):
        self.logger = logging.getLogger('iterativeWGCNA.Network')
        self.genes = genes
        self.classifiedGenes = self.genes.get_classified_genes()
        self.profiles = self.genes.profiles
        # self.graph = None
        self.args = args
        self.eigengenes = eigengenes
        self.modules = self.genes.get_modules()
        self.moduleColors = {}
        self.generate_module_colors()
        self.genes.assign_colors(self.moduleColors)

        self.weightMatrix = None # membership weight matrix
        # self.generate_weight_matrix()

        self.adjacency = None
        self.weightedAdjacency = None # weighted by shared membership
        # self.generate_weighted_adjacency()


    def summarize(self):
        '''
        generate summary figs and data
        '''
        self.plot_eigengene_network()
        if self.args.generateNetworkSummary is not None:
            self.__plot_summary_views()


    def generate_module_colors(self):
        '''
        generate module color map
        '''
        for m in self.modules:
            # TODO probably a better way to do this, but...
            # since we may need > 50 colors, random it is
            self.moduleColors[m] = '#' + '%06X' % randint(0, 0xFFFFFF)
        self.moduleColors["UNCLASSIFIED"] = '#D3D3D3'


    def plot_eigengene_network(self):
        '''
        wrapper for plotting the eigengene network to pdf
        '''
        grdevices().pdf("eigengene-network.pdf")
        manager = WgcnaManager(self.eigengenes.matrix, self.args.wgcnaParameters)
        manager.plot_eigengene_network()
        grdevices().dev_off()


    def plot_network_summary(self, genes, title, filename):
        '''
        wrapper for WGCNA summary network view (heatmap + dendrogram)
        '''

        expression = self.profiles.gene_expression(genes)
        colors = self.genes.get_module_colors(genes)
        manager = WgcnaManager(expression, self.args.wgcnaParameters)
        grdevices().pdf(filename)
        manager.plot_network_overview(colors, title)
        grdevices().dev_off()


    def __plot_summary_views(self):
        '''
        plot summary heatmaps/views
        '''
        if self.args.generateNetworkSummary == 'input' \
           or self.args.generateNetworkSummary == 'all':
            self.plot_network_summary(self.genes.get_genes(),
                                      "All Genes (incl. unclassified)",
                                      "input-overview.pdf")

        if self.args.generateNetworkSummary == 'network' \
           or self.args.generateNetworkSummary == 'all':
            self.plot_network_summary(self.classifiedGenes,
                                      "Network (classified genes)",
                                      "network-overview.pdf")


    def generate_weighted_adjacency(self):
        '''
        gene x gene weight matrix with matrix[r][c] = 1
        if genes r & c are in the same module; for
        weighted graph viz and to simply
        calc of in/out degree
        '''

        manager = WgcnaManager(self.profiles.gene_expression(self.classifiedGenes),
                               self.args)
        manager.adjacency(True, True, True) # signed, but filter negatives & self-refs
        self.adjacency = manager.adjacencyMatrix
        self.weightedAdjacency = self.adjacency

        for m in self.modules:
            members = self.genes.get_module_members(m)
            for r in members:
                indexR = self.weightedAdjacency.names.index(r)
                self.logger.debug(r)
                self.logger.debug(self.weightedAdjacency.names)
                self.logger.debug("gene: " + r + " index " + str(indexR))

                for c in members:
                    adj = self.adjacency.rx(r, c)
                    self.logger.debug(adj)

                    if adj > 0:
                        self.logger.debug(self.weightedAdjacency.names)
                        indexC = self.weightedAdjacency.names.index(c)
                        self.weightedAdjacency[indexR][indexC] = adj + 1
                        self.logger.debug(self.weightedAdjacency.rx(r, c))

                    sys.exit(1)
