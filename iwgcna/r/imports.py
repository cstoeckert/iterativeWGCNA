'''
imports from R; wrapped in functions
to ensure warning messages go to the R log
'''
from rpy2.robjects.packages import importr, SignatureTranslatedAnonymousPackage
from .snippets import UTIL_FUNCTIONS

r_utils = SignatureTranslatedAnonymousPackage(UTIL_FUNCTIONS, 'r_utils')

def base():
    return importr('base')


def wgcna():
    return importr('WGCNA')


def igraph():
    return importr('igraph')


def stats():
    return importr('stats')


def graphics():
    return importr('graphics')


def grdevices():
    return importr('grDevices')



