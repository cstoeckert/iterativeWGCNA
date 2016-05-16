'''
manage eigengene matrices
'''

import rpy2.robjects as ro
from ..utils.io import write_data_frame
from ..r.imports import base, stats

def write(ematrix, isFinal):
    '''
    writes the eigengene matrix to file
    '''
    fileName = 'eigengenes-final.txt' if isFinal else 'eigengenes.txt'   
    write_data_frame(ematrix, fileName, 'Module')


def equal(e1, e2):
    '''
    check if 2 eigengenes are "equal" (correlation ~ 1)
    '''
    correlation = base().as_data_frame(stats().cor(base().t(e1), base().t(e2)))
    correlation = round(correlation.rx(1, 1)[0], 1)
    return correlation == 1


def load_from_file(fileName):
    '''
    loads eigengenes from file into an R DataFrame
    '''
    return ro.DataFrame.from_csvfile(fileName, sep='\t', header=True, row_names=1)


def extract_modules(ematrix, modules):
    return ematrix.rx(ro.StrVector(modules), True)


