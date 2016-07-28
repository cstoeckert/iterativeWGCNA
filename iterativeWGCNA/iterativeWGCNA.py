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
        self.profiles = exprData
        self.args = args
        self.profiles = None

        # create working directory
        io.create_dir(self.args.workingDir)

        self.__initialize_R()

        # initalize logger
        logger = log.initialize(args.workingDir)

        log.parameters(args)
        if args.verbose:
            log.parameters(args)


    def load_expression_profiles(self):
          # gives a weird R error that I'm having trouble catching
          # when it fails
        # TODO: identify the exact exception 
        try:
            data = io.read_data(args.inputFile)
        except:
            logger.error("Unable to open input file: " + args.inputFile)
            sys.exit(1)

        log.input_data(args.inputFile, data.ncol, data.nrow)
        
            
        
    def __initialize_R(self):
        '''
        initialize R workspace and logs
        '''
        
        initialize_r_workspace(self.args.workingDir, self.args.allowWGCNAThreads)
