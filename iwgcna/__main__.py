#!/usr/bin/env python2.7
"""
Perform iterative WGCNA analysis

python dependencies:
  * rpy2

R dependencies:
  * igraph
  * WGCNA
  * lattice
"""

from . import parse_command_line_args
from .iterativeWGCNA import main

if __name__ == '__main__':
    CML_ARGS = parse_command_line_args()
    main(CML_ARGS)

__author__ = 'Emily Greenfest-Allen'
__copyright__ = 'Copyright 2016, University of Pennsylvania'
