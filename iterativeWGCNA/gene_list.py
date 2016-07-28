# pylint: disable=invalid-name

'''
manage gene list
'''

import logging
from collections import OrderedDict

# from .expression import Expression
from .gene_properties import GeneProperties

class Genes(object):
    '''
    track input genes and their properties, including
    expression profiles, module membership, eigengene
    connectivity
    '''

    def __init__(self, exprData):
        '''
        initialize an OrderedDict of genes
        from the row.names of the expression
        data set
        '''
        genes = exprData.genes()

        self.genes = OrderedDict((geneId, None) for geneId in genes)
        self.size = len(self.genes)

        for geneId in genes:
            props = GeneProperties(exprData.gene_expression(geneId))
            self.genes[geneId] = props
