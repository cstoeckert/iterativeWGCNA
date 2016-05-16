'''
calculate and manage eigengene connectivity lists (kME)
'''

from collections import OrderedDict
import rpy2.robjects as ro
from .expression import get_member_expression
from .r.imports import wgcna, stats, base
from .io.utils import write_data_frame

def initialize(data):
    '''
    initialized eigengene connectivity (kME)
    dictionary
    gene list comes from input data row names (DATA.rownames)
    all kMEs are initially NaN
    an ordered dictionary is used to keep values in the same
    order as input data
    '''
    kME = OrderedDict((gene, float('NaN')) for gene in data.rownames)
    return kME


def calculate(expr, eigengene, calculateP):
    '''
    calculates eigengene connectivity kme
    between an eigengene and expression data set
    '''
    if calculateP:
        correlation = wgcna().corAndPvalue(base().t(expr), base().t(eigengene))
    else:
        correlation = base().as_data_frame(stats().cor(base().t(expr), base().t(eigengene)))
    return correlation


def update(kME, data, membership, eigengenes):
    '''
    finds new module membership and updates eigengene
    connectivity (kME)
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
