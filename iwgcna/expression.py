'''
functions for manipulating expression matrices
'''

import rpy2.robjects as ro
from .r.imports import r_utils

def get_member_expression(module, expr, membership):
    '''
    subsets expression data
    returning expression for only members of the specified module
    '''
    return r_utils.extractMembers(module, expr, ro.ListVector(membership))


def get_residuals(expr, membership):
    '''
    subsets expression data
    returning expression for only residuals to the fit
    '''
    if membership is None:
        return expr
    else:
        return get_member_expression('UNCLASSIFIED', expr, membership)


def remove_residuals(expr, membership):
    '''
    subsets expression data
    removing residuals to the fit
    '''
    if membership is None:
        return expr
    else:
        return r_utils.removeUnclassified(expr, ro.ListVector(membership))
