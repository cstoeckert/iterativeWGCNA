#!/usr/bin/env python2.7
"""
Perform iterative WGCNA analysis

python dependencies:
  * rpy2
  * matplotlib

R dependencies:
  * WGCNA
"""

from .cmlargs import parse_command_line_args
from .iterativeWGCNA import IterativeWGCNA

if __name__ == '__main__':
    args = parse_command_line_args()
    alg = IterativeWGCNA(args)
    alg.run()

__author__ = 'Emily Greenfest-Allen'
__copyright__ = 'Copyright 2016, University of Pennsylvania'
