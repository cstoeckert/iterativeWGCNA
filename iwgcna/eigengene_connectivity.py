'''
calculate and manage eigengene connectivity lists (kME)
'''

from collections import OrderedDict
import rpy2.robjects as ro

from .expression import get_member_expression

from .r.imports import wgcna, stats, base
from .io.utils import write_data_frame

class EigeneneConnectivityList(object):
    ''' 
    track, calculate, and manage eigengene 
    connectivity (Kme) for sequence features
    '''

    def __init__(self, features):
        '''
        initialize an OrderedDict with one entry per feature
        all features have initial Kme of NaN
        '''
        self.values = OrderedDict((f, float('NaN')) for f in features)
        self.size = len(features)
        self.iteration = None
        return None

    def calculate_kme(self, expr, eigengene, calculateP):
        '''
        calculates eigengene connectivity
        between an eigengene and expression data set
        '''
        if calculateP:
            correlation = wgcna().corAndPvalue(base().t(expr), base().t(eigengene))
        else:
            correlation = base().as_data_frame(stats().cor(base().t(expr), base().t(eigengene)))
            return correlation


def update(kME, data, membership, eigengenes):
    '''
    updates eigengene connectivity (kME)
    for each module, extracts the member subset from the
    expression data and calculates the kME between the module
    eigengene and each member
    '''

    for module in eigengenes.rownames:
        moduleEigengene = eigengenes.rx(module, True)
        moduleMemberExpression = get_member_expression(module, data, membership)
        memberKME = calculate(moduleMemberExpression, moduleEigengene, False)
        for gene in memberKME.rownames:
            kME[gene] = round(memberKME.rx(gene, 1)[0], 2)

    return kME


def write(iteration, kME):
    '''
    writes the eigengene connectivity (kME)
    dictionary to file
    '''
    df = ro.DataFrame(kME)
    df.rownames = (iteration)
    write_data_frame(df, 'eigengene-connectivity.txt', 'Iteration')
