# pylint: disable=invalid-name
# pylint: disable=no-self-use
# pylint: disable=too-many-instance-attributes
# pylint: disable=redefined-variable-type

'''
manage network (for summaries)
'''

from __future__ import print_function
from __future__ import with_statement

import logging

from collections import OrderedDict

import rpy2.robjects as ro

from .colors import Colors
from .wgcna import WgcnaManager
from .r.manager import RManager
from .r.imports import grdevices, base, rsnippets

from .eigengenes import Eigengenes

class Network(object):
    '''
    a network contains the classified genes
    and their assignments
    hash of modules and properties
    for summary purposes
    '''

    def __init__(self, args):
        self.logger = logging.getLogger('iterativeWGCNA.Network')
        self.args = args
        self.eigengenes = None

        self.genes = None
        self.modules = None
        self.classifiedGenes = None
        self.profiles = None
        self.kME = None
        self.membership = None

        # self.graph = None
        self.geneColors = None
        self.adjacency = None
        self.weightedAdjacency = None # weighted by shared membership


    def build(self, genes, eigengenes):
        '''
        build from iterativeWGCNA result in memory
        '''
        self.eigengenes = eigengenes
        self.genes = genes.get_genes()
        self.classifiedGenes = genes.get_classified_genes()
        self.profiles = genes.profiles
        self.kME = genes.get_gene_kME() # TODO -- fix this -- this function has changed
        self.membership = genes.get_gene_membership()

        self.modules = genes.get_modules()
        self.__initialize_module_properties()

        self.__assign_colors()
        self.__generate_weighted_adjacency()


    def __assign_colors(self):
        self.__generate_module_colors()
        self.geneColors = OrderedDict((gene, None) for gene in self.genes)
        self.assign_gene_colors()


    def __initialize_module_properties(self):
        '''
        transform module list into dict with placeholders
        for color, kIn, and kOut
        and values for size
        '''
        self.modules = OrderedDict((module,
                                    {'color':None,
                                     'kIn':0,
                                     'kOut':0,
                                     'size':self.__get_module_size(module),
                                     'density':0}) \
                                       for module in self.modules)
        self.modules.update({'UNCLASSIFIED': {'color':None, 'kIn': 0,
                                              'kOut': 0,
                                              'size': self.__get_module_size('UNCLASSIFIED'),
                                              'density': 0.0}})


    def build_from_file(self, profiles, adjacency=True):
        '''
        initialize Network from iterativeWGCNA output found in path
        '''
        self.profiles = profiles
        self.genes = self.profiles.genes()

        # when membership is loaded from file, modules and classified
        # genes are determined as well
        self.__load_membership_from_file(self.args.preMerge)
        self.__load_kme_from_file(self.args.preMerge)
        warning("done")
        
        self.eigengenes = Eigengenes()
        self.eigengenes.load_matrix_from_file("eigengenes-final.txt")
        if adjacency:
            self.__initialize_module_properties()
            self.__assign_colors()
            self.__generate_weighted_adjacency()


    def __load_membership_from_file(self, preMerge):
        '''
        loads membership assignments from file
        and assembles list of classified genes
        and determines list of unique modules
        '''
        fileName = "membership.txt"
        membership = ro.DataFrame.from_csvfile(fileName, sep='\t',
                                               header=True, row_names=1, as_is=True)

        finalIndex = membership.names.index('final')
        if preMerge:
            finalIndex = finalIndex - 1 

        self.membership = OrderedDict((gene, None) for gene in self.genes)
        self.classifiedGenes = []
        self.modules = []

        for g in self.genes:
            module = membership.rx(g, finalIndex)[0]
            self.modules.append(module)
            self.membership[g] = module
            if module != 'UNCLASSIFIED':
                if 'p' not in module:
                    self.logger.debug("Module: " + module + "; Gene: " + g)
                self.classifiedGenes.append(g)

        self.modules = list(set(self.modules)) # gets unique list of modules


    def __load_kme_from_file(self, preMerge):
        '''
        loads kME to assigned module from file
        '''
        fileName = "eigengene-connectivity.txt"
        kME = ro.DataFrame.from_csvfile(fileName, sep='\t',
                                        header=True, row_names=1, as_is=True)

        finalIndex = kME.names.index('final')
        if preMerge:
            finalIndex = finalIndex - 1
            
        self.kME = OrderedDict((gene, None) for gene in self.genes)
        for g in self.genes:
            self.kME[g] = kME.rx(g, finalIndex)[0]


    def summarize_network(self):
        '''
        generate summary figs and data
        '''
        self.__plot_eigengene_overview()
        self.__plot_summary_views()
        self.__summarize_network_modularity()
        self.__write_module_summary()


    def summarize_module(self, module):
        '''
        generate summare info for the specified module
        '''
        self.__plot_module_overview(module)


    def __plot_module_overview(self, module):
        '''
        plot heatmap, dendrogram, and eigengene
        '''
        grdevices().pdf(module + "-summary.pdf")
        self.plot_module_eigengene(module)
        self.plot_module_kME(module)
        self.plot_module_heatmap(module)
        grdevices().dev_off()


    def plot_module_heatmap(self, module):
        '''
        plot module heatmap
        '''

        members = self.__get_module_members(module)
        expression = self.profiles.gene_expression(members)
        manager = RManager(expression, None)
        manager.heatmap()


    def plot_module_eigengene(self, module):
        '''
        barchart illustrating module eigengene
        '''
        eigengene = self.eigengenes.get_module_eigengene(module)

        params = {}
        params['height'] = base().as_numeric(eigengene)

        limit = max(abs(base().max(eigengene)[0]), abs(base().min(eigengene)[0]))
        ylim = [-1 * limit, limit]
        params['ylim'] = ro.IntVector(ylim)

        colors = ["red" if e[0] > 0 else "blue" for e in eigengene]
        params['col'] = ro.StrVector(colors)

        params['border'] = ro.NA_Logical
        params['las'] = 2
        params['names.arg'] = ro.StrVector(self.eigengenes.samples())
        params['cex.names'] = 0.6
        params['main'] = "Eigengene: " + module
        manager = RManager(eigengene, params)
        manager.barchart()


    def plot_module_kME(self, module):
        '''
        plots module eigengene connectivity (kME)
        '''
        members = self.__get_module_members(module)
        kME = [kME for gene, kME in self.kME.items() if gene in members]

        manager = RManager(kME, None)
        manager.histogram(self.args.wgcnaParameters['minKMEtoStay'],
                          {'main':"Member kME: " + module,
                           'xlab': "Eigengene Connectivity (kME)",
                           'ylab': "N Genes",
                           'breaks': base().seq(0, 1, 0.1)})



    def __generate_module_colors(self):
        '''
        generate module color map
        using standard colors for first
        few modules then random for anything more than first 25
        '''

        colors = Colors()
        n = 0 # counter for module num
        for m in self.modules:
            n = n + 1
            color = colors.assign_color(n)
            self.modules[m]['color'] = color

        self.modules["UNCLASSIFIED"]['color'] = '#D3D3D3'


    def assign_gene_colors(self):
        '''
        assign colors to genes according to module membership
        '''
        for g in self.genes:
            self.geneColors[g] = self.modules[self.membership[g]]['color']


    def get_gene_colors(self, targetGenes):
        '''
        retrieve colors for specified gene list
        '''
        colors = OrderedDict((gene, color) for gene, color in self.geneColors.items() \
                                 if gene in targetGenes)
        return colors


    def get_gene_membership(self, targetGenes):
        '''
        retrieve membership for specified gene list
        '''
        colors = OrderedDict((gene, membership) for gene, membership in self.membership.items() \
                                 if gene in targetGenes)
        return colors


    def __plot_eigengene_overview(self):
        '''
        plots eigengene graphs to single pdf
        '''
        grdevices().pdf("eigengene-overview.pdf")
        self.plot_eigengene_network()
        self.plot_eigengene_heatmap()
        grdevices().dev_off()


    def plot_eigengene_network(self):
        '''
        wrapper for plotting the eigengene network to pdf
        '''
        manager = WgcnaManager(self.eigengenes.matrix, None)
        manager.plot_eigengene_network()


    def plot_eigengene_heatmap(self):
        '''
        plot a heatmap of the eigengenes
        '''
        manager = RManager(self.eigengenes.matrix, None)
        # manager.heatmap(params={'scale':'none'})
        manager.heatmap()

        
    def plot_network_summary(self, genes, title, filename):
        '''
        wrapper for WGCNA summary network view (heatmap + dendrogram)
        '''

        expression = self.profiles.gene_expression(genes)
        membership = self.get_gene_membership(genes)
        # colors = self.get_gene_colors(genes)
        
        manager = WgcnaManager(expression, self.args.wgcnaParameters)
        manager.set_module_colors(self.modules)
        
        grdevices().pdf(filename)
        manager.plot_network_heatmap(membership, title) # plot_network_overview(colors, title)
        grdevices().dev_off()


    def __plot_summary_views(self):
        '''
        plot summary heatmaps/views
        '''
        if self.args.generateNetworkSummary == 'input' \
           or self.args.generateNetworkSummary == 'all':
            self.plot_network_summary(self.genes,
                                      "All Genes (incl. unclassified)",
                                      "input-block-diagram.pdf")

        if self.args.generateNetworkSummary == 'network' \
           or self.args.generateNetworkSummary == 'all':
            self.plot_network_summary(self.classifiedGenes,
                                      "Network (classified genes)",
                                      "network-block-diagram.pdf")


    def __get_module_members(self, targetModule):
        '''
        get genes in targetModule
        '''
        return [gene for gene, module in self.membership.items() if module == targetModule]


    def __generate_weighted_adjacency(self):
        '''
        gene x gene weight matrix with matrix[r][c] = 1
        if genes r & c are in the same module; for
        weighted graph viz and to simply
        calc of in/out degree
        '''

        manager = WgcnaManager(self.profiles.gene_expression(self.classifiedGenes),
                               self.args.wgcnaParameters)
        manager.adjacency('signed', True, True) # signed, but filter negatives & self-refs
        self.adjacency = base().as_data_frame(manager.adjacencyMatrix)
        self.weightedAdjacency = self.adjacency

        for m in self.modules:
            if m == 'UNCLASSIFIED':
                continue

            members = self.__get_module_members(m)
            for r in members:
                indexR = self.weightedAdjacency.names.index(r)
                for c in members:
                    if r == c:
                        continue
                    indexC = self.weightedAdjacency.names.index(c)
                    adj = self.adjacency[indexR][indexC]
                    if adj > 0:
                        self.weightedAdjacency[indexR][indexC] = adj + 1


    def calculate_degree_modularity(self, targetModule):
        '''
        calculates in degree (kIn) and out degree (kOut)
        for the target module
        '''
        members = self.__get_module_members(targetModule)

        degree = rsnippets.degree(self.adjacency, ro.StrVector(members),
                                  self.args.edgeWeight)
        self.modules[targetModule]['kIn'] = int(degree.rx2('kIn')[0])
        self.modules[targetModule]['kOut'] = int(degree.rx2('kOut')[0])
        size = self.modules[targetModule]['size']
        self.modules[targetModule]['density'] = float(self.modules[targetModule]['kIn'])/(float(size) * (float(size) - 1.0)/2.0)


    def __summarize_network_modularity(self):
        '''
        summarize modularity of network
        calculates in degree (kIn) and out degree (kOut) per module
        '''
        for m in self.modules:
            if m != 'UNCLASSIFIED':
                self.calculate_degree_modularity(m)


    def __get_module_size(self, targetModule):
        '''
        return # of module members for target module
        '''
        return len(self.__get_module_members(targetModule))


    def __write_module_summary(self):
        '''
        writes network modularity to a file
        '''
        fileName = 'module-summary.txt'
        with open(fileName, 'w') as f:
            header = ('module', 'size', 'color',
                      'kOut', 'avg_node_kOut',
                      'kIn',
                      'module_density', 'kIn2kOut_ratio')
            print('\t'.join(header), file=f)
            for m in self.modules:
                kIn = self.modules[m]['kIn']
                kOut = self.modules[m]['kOut']
                size = self.modules[m]['size']
                avg_kOut = "{0:0.0f}".format(round(float(kOut) / float(size)))
                ratio = "{0:0.2f}".format(float(kIn) / float(kOut)) \
                         if kOut != 0 else 'NA'

                print(m,
                      size,
                      self.modules[m]['color'],
                      kOut, avg_kOut,
                      kIn,
                      "{0:0.1f}".format(self.modules[m]['density']),
                      ratio,
                      sep='\t', file=f)


    def export_cytoscape_json(self):
        '''
        creates and saves a cytoscape json file
        '''
        filterLevel = self.args.edgeWeight

        if self.weightedAdjacency is None:
            self.__generate_weighted_adjacency()

