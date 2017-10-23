# pylint: disable=invalid-name
# pylint: disable=bare-except
# pylint: disable=broad-except
# pylint: disable=too-many-instance-attributes
'''
plottter: executes plotting/handles logging and workspace
'''

from __future__ import print_function

import logging
import sys
import os
import argparse
from time import strftime
from os import getcwd

import rpy2.robjects as ro
from .io.utils import create_dir, read_data, warning, write_data_frame
from .r.manager import RManager
from .r.imports import base, circlize, reshape2, colorRamps


class PlotManager(RManager):
    '''
    manages plotting
    '''
    def __init__(self, logPrefix='plot'):
        self.__initialize_log(logPrefix)
        self.logger.info(strftime("%c"))
        RManager.__init__(self, None, None)
        
        self.args = self.__parse_command_line_args()
        self.__initialize_R(logPrefix)


    def __initialize_log(self, prefix='plot'):
        '''
        initialize log by setting path and file format
        '''
        logging.basicConfig(filename=self.args.workingDir + '/' + prefix + '.log',
                            filemode='w', format='%(levelname)s: %(message)s',
                            level=logging.DEBUG)

        logging.captureWarnings(True)
        self.logger = logging.getLogger(__name__)


    def __initialize_R(self, prefix='plot'):
        '''
        initialize R workspace and logs
        '''
        # set working directory
        base().setwd(self.args.dir)

        # suppress warnings
        ro.r['options'](warn=-1)

        # r log
        logFile = prefix + '.log'
        rLogger = base().file(logFile, open='wt')
        base().sink(rLogger, type=base().c('output', 'message'))


        
    def __parse_command_line_args(self):
        '''
        parse command line args for plots
        '''
        
        parser = argparse.ArgumentParser(prog='plot iterativeWGCNA results',
                                     description="generate graphical results from interative WGCNA analysis",
                                     formatter_class=argparse.RawTextHelpFormatter)

        parser.add_argument('--dir',
                            help="directory containing the iterativeWGCNA output to be plotted",
                            metavar='<output dir>',
                                default=getcwd())
        
        parser.add_argument('--pass',
                            help="(optional) pass number (default = 1)",
                            default=1,
                            type=int)

        parser.add_argument('--iteration',
                            help="(optional) iteration number (default = 1); together with --pass specify iteration against which final result is to be compared; default is the first iteration (pass1/iteration1)",
                            default=1,
                            type=int)

        parser.add_argument('-v', '--verbose',
                            help="print status messages",
                            action='store_true')
        
        return parser.parse_args()
    
    def load_blockwiseWGCNA_result(self):
        filePath = os.path.join(self.args.dir, 'pass' + str(self.args.pass),
                                'i' + str(self.args.iteration))
        
