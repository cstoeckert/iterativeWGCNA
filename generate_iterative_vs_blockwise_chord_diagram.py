#!/usr/bin/env python

'''
Use the circlize R package to plot a chord diagram comparing clustering from blockwise WGCNA (result from pass1, iteration 1) to final iterativeWGCNA clustering
'''

from iterativeWGCNA.cmlargs import parse_plotting_command_line_args
from iterativeWGCNA.plot import PlotManager


if __name__ == '__main__':
    manager = PlotManager('wgcna_comparision_chord_diagram')
    manager.load_blockwiseWGCNA_result()
    manager.load_iterativeWGCNA_result()
    manager.generate_ribbons()
    manager.write_ribbons()
    manager.chord_diagram()

__author__ = 'Emily Greenfest-Allen based on original script by JP Cartailler'
__copyright__ = 'Copyright 2016, University of Pennsylvania'

