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

import logging
from time import strftime

from iterativeWGCNA.cmlargs import parse_summary_command_line_args
from iterativeWGCNA.iterativeWGCNA import IterativeWGCNA
from iterativeWGCNA.network import Network

def add_wgcna_params(args):
    '''
    network functions expect a 'wgcnaParameters' item
    in the argparse Namespace, so adding it in
    '''
    args.wgcnaParameters = {}
    args.wgcnaParameters['power'] = args.power
    args.wgcnaParameters['minKMEtoStay'] = args.minKMEtoStay
    args.wgcnaParameters['type'] = 'signed' if args.signed \
                                   else 'unsigned'
    return args


if __name__ == '__main__':
    cmlArgs = parse_summary_command_line_args()
    alg = IterativeWGCNA(cmlArgs, True) # initialize workspace, logs and load profiles
    logger = logging.getLogger('iterativeWGCNA.SummarizeResults')

    try:
        cmlArgs = add_wgcna_params(cmlArgs)
        network = Network(cmlArgs)
        network.build_from_file(alg.profiles)

        if cmlArgs.generateNetworkSummary:
            network.summarize_network()

        if cmlArgs.summarizeModule is not None:
            network.summarize_module(cmlArgs.summarizeModule)

        logger.info('Summary: SUCCESS')
    except Exception:
        if logger is not None:
            logger.exception('Summary: FAIL')
        else:
            raise
    finally:
        if logger is not None:
            logger.info(strftime("%c"))

__author__ = 'Emily Greenfest-Allen'
__copyright__ = 'Copyright 2016, University of Pennsylvania'
