#!/usr/bin/env python

'''Rerun final module merge'''

# Installation workaround - see README
# import readline
#pylint: disable=invalid-name

from iterativeWGCNA.iterativeWGCNA import IterativeWGCNA
from iterativeWGCNA.cmlargs import parse_command_line_args

if __name__ == '__main__':
    args = parse_command_line_args(program='iterativeWGCNA: Adjust Merge',
                                   description='recompute final module merge from existing output')
    alg = IterativeWGCNA(args, report="merge")
    alg.merge_close_modules_from_output()

__author__ = 'Emily Greenfest-Allen'
__copyright__ = 'Copyright 2018, University of Pennsylvania'


