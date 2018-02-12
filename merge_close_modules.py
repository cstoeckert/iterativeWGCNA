#!/usr/bin/env python

'''Rerun final module merge'''

# Installation workaround - see README
# import readline
#pylint: disable=invalid-name

import argparse
from os import getcwd

from iterativeWGCNA.iterativeWGCNA import IterativeWGCNA
from iterativeWGCNA.cmlargs import restricted_float


if __name__ == '__main__':
    args = parse_command_line_args()
    alg = IterativeWGCNA(args, report="merge")
    alg.merge_close_modules_from_output()

__author__ = 'Emily Greenfest-Allen'
__copyright__ = 'Copyright 2018, University of Pennsylvania'



def parse_command_line_args():
    '''
    parses command line args
    '''

    parser = argparse.ArgumentParser(prog='iterativeWGCNA-Adjust Merge',
                                     description="recomputes final module merge, given an input file and path to iterativeWGCNA output")

    parser.add_argument('-i', '--inputFile',
                        metavar='<gene expression file>',
                        help="full path to input gene expression file; "
                        + "if full path is not provided,\n"
                        + "assumes the file is in the working directory\n;"
                        + "see below for input file format",
                        required=True)

    parser.add_argument('-o', '--workingDir',
                        help="R working directory; where output will be saved",
                        metavar='<output dir>', 
                        default=getcwd())

    parser.add_argument('-d', '--debug',
                        help="print additional messages for debugging",
                        action='store_true')
    
    parser.add_argument('-v', '--verbose',
                        help="print status messages",
                        action='store_true')

    parser.add_argument('-f', '--finalMergeCutHeight',
                        help="cut height for final merge (after iterations are assembled)",
                        default=0.05,
                        metavar='<cut height>',
                        type=restricted_float)

    return parser.parse_args()
