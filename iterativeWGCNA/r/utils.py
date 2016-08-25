'''
initialize or manage the R workspace
'''

from .imports import base, wgcna, grdevices
import rpy2.robjects as ro

def initialize_r_workspace(workingDir, allowThreads):
    '''
    initialize the r working environment
    and r log
    also allow WGCNA threading if specified in
    command line arguments
    '''

    # set working directory
    base().setwd(workingDir)

    # suppress warnings
    ro.r['options'](warn=-1)

    # r log
    rLogger = base().file('R.log', open='wt')
    base().sink(rLogger, type=base().c('output', 'message'))

    if allowThreads:
        wgcna().enableWGCNAThreads()


def barplot(valueMap, fileName):
    '''
    takes a dict of label=value pairs and generates
    a bar plot
    '''
    grdevices().png(file=fileName, width=512, height=512)
    categories = list(valueMap.keys())
    values = list(valueMap.values())
    grdevices().dev_off()
