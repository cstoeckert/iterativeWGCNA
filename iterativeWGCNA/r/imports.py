'''
imports from R; wrapped in functions
to ensure warning messages go to the R log
'''
from rpy2.robjects.packages import importr, SignatureTranslatedAnonymousPackage
from .snippets import FUNCTIONS

rsnippets = SignatureTranslatedAnonymousPackage(FUNCTIONS, 'rsnippets')

def base():
    return importr('base')


def wgcna():
    return importr('WGCNA')


def stats():
    return importr('stats')


def graphics():
    return importr('graphics')


def grdevices():
    return importr('grDevices')


def pheatmap():
    return importr('pheatmap')
