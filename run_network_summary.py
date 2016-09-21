#!/usr/bin/env python2.7

# pylint: disable=invalid-name

'''
Convenience wrapper for summarizing output
from a run of iterativeWGCNA
can either summarize the final network or
specify an module or specific iteration
'''

# Installation workaround - see README
# import readline

from iterativeWGCNA.cmlargs import parse_summary_command_line_args
from iterativeWGCNA.iterativeWGCNA import IterativeWGCNA
from iterativeWGCNA.network import Network

if __name__ == '__main__':
    args = parse_summary_command_line_args()
    alg = IterativeWGCNA(args) # initialize workspace, logs and load profiles
    network = Network(args)
    network.build_from_file(alg.profiles)
    network.summarize()

__author__ = 'Emily Greenfest-Allen'
__copyright__ = 'Copyright 2016, University of Pennsylvania'
