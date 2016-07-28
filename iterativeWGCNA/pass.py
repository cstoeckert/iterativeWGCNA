# pylint: disable=invalid-name
# pylint: disable=unused-import
'''
a pass of iterativeWGCNA
'''

import logging

from .genes import Genes
from .iteration import Iteration
from .io.utils import warning
from .r.utils import initialize_r_workspace

#from .analysis import calculate_kME
#from .r.imports import wgcna, stats, base, rsnippets

class Pass(object):
    '''
	a pass of iterativeWGCNA: iterative application of WGCNA
    until no more genes are prunned from the network
    '''

    def __init__(self, exprData, args):
        self.logger = logging.getLogger('iterativeWGCNA.Pass')

        self.genes = Genes(exprData)
        self.modules = {}
        self.iteration = None
        self.args = args
    
