# pylint: disable=invalid-name
# pylint: disable=unused-import
'''
main application
'''

import logging
from collections import OrderedDict
#from collections import Counter

import rpy2.robjects as ro

from .genes import Genes

#from .analysis import calculate_kME
#from .r.imports import wgcna, stats, base, rsnippets
#from .io.utils import write_data_frame

class iterativeWGCNA(object):
    '''
    main application
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
        initialize_r_workspace(args.workingDir, args.allowWGCNAThreads)

        logger = log.initialize(args.workingDir)
        log.parameters(args)
