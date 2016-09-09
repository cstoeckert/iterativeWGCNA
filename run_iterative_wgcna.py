#!/usr/bin/env python2.7

'''Convenience wrapper for running iterativeWGCNA directly from source tree.'''

# Installation workaround - see README
# import readline

from iterativeWGCNA.cmlargs import parse_command_line_args
from iterativeWGCNA.iterativeWGCNA import IterativeWGCNA

if __name__ == '__main__':
    args = parse_command_line_args()
    alg = IterativeWGCNA(args)
    alg.run()

__author__ = 'Emily Greenfest-Allen'
__copyright__ = 'Copyright 2016, University of Pennsylvania'
