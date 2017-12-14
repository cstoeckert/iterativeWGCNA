'''
functions in support of data analysis
'''

# TODO move to RManager or wgcnaManager

from .r.imports import wgcna, stats, base

def calculate_kME(expr, eigengene, calculateP):
    '''
    calculates eigengene connectivity
    between an eigengene and expression data set
    '''
    if calculateP:
        correlation = wgcna().corAndPvalue(base().t(expr), base().t(eigengene))
    else:
        correlation = base().as_data_frame(stats().cor(base().t(expr), base().t(eigengene)))
    return correlation
