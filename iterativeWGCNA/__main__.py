#!/usr/bin/env python2.7
"""
Perform iterative WGCNA analysis

python dependencies:
  * rpy2
  * matplotlib

R dependencies:
  * igraph
  * WGCNA
"""

from . import parse_command_line_args
from .iterativeWGCNA import iterativeWGCNA

if __name__ == '__main__':
    iterativeWGCNA(parse_command_line_args())

__author__ = 'Emily Greenfest-Allen'
__copyright__ = 'Copyright 2016, University of Pennsylvania'
