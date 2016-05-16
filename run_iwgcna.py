#!/usr/bin/env python2.7

'''Convenience wrapper for running wgcna directly from source tree.'''

from iwgcna.iterativeWGCNA import main
from iwgcna import parse_command_line_args

if __name__ == '__main__':
    CML_ARGS = parse_command_line_args()
    main(CML_ARGS)
