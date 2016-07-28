# pylint: disable=invalid-name
# pylint: disable=unused-import
'''
result from a single iteration
'''

import logging

from .genes import Genes
from .io.utils import warning
import .r.utils as rutils

#from .analysis import calculate_kME
#from .r.imports import wgcna, stats, base, rsnippets

class Iteration(object):
    '''
	an iteration of iterativeWGCNA, including both
	data and 
    '''

    def __init__(self, exprData, args):
        self.logger = logging.getLogger('iterativeWGCNA')
        # genes
        self.genes = Genes(exprData)
        self.modules = {}
        self.iteration = None
        self.args = args

        # initialize R workspace and logs
        io.create_dir(args.workingDir)
        rutils.initialize_r_workspace(args.workingDir, args.allowWGCNAThreads)

        logger = log.initialize(args.workingDir)
        log.parameters(args)
