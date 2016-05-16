'''
manage eigengene matrices
'''

import rpy2.robjects as ro

from .io.utils import write_data_frame
from .r.imports import base, stats

def write(ematrix, isFinal):
    '''
    writes the eigengene matrix to file
    '''
    fileName = 'eigengenes-final.txt' if isFinal else 'eigengenes.txt'
    write_data_frame(ematrix, fileName, 'Module')

def similarity(ematrix, module):
    sim = base().as_data_frame(stats().cor(base().t(ematrix), \
             base().t(ematrix.rx(module, True))))    
    return sim
    
def correlation(e1, e2):
    '''
    calculate correlation between two eigengenes
    '''
    cor = base().as_data_frame(stats().cor(base().t(e1), base().t(e2)))
    cor = round(cor.rx(1, 1)[0], 1)
    return cor


def equal(e1, e2, threshold=0.0):
    '''
    check if 2 eigengenes are "equal" (1 - correlation <= threshold)
    '''
    cor = correlation(e1, e2)
    return 1.0 - cor <= threshold


def load_from_file(fileName):
    '''
    loads eigengenes from file into an R DataFrame
    '''
    return ro.DataFrame.from_csvfile(fileName, sep='\t', header=True, row_names=1)


def extract_modules(ematrix, modules):
    return ematrix.rx(ro.StrVector(modules), True)


