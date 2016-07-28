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
    connectivity for sequence features
    '''

    def __init__(self, features):
        '''
        initialize an OrderedDict with one entry per feature
        all features have initial kme of NaN
        '''
        self.values = OrderedDict((f, float('NaN')) for f in features)
        self.size = len(features)
        self.iteration = None
        return None

    def __update(self, feature, kME):
        '''
        update value of 'feature' to kme
        do not add new featuers
        '''
        if feature in self.values:
            self.values = kME
            return True
        else:
            return False
        

    # TODO def calculate(expr, eigengene, calculateP):
    def __calculate(epression, eigengene, calculateP):
        
    
    def update(self, module):
        '''
        updates feature eigengene connectivity by calculating 
        the kME between the module eigengene and member 
        expression profiles

        :param module: a module object, containing the features, eigengene, etc.
        '''

        memberKME = self.__calculate(module.expression, module.eigengene, False)
        for feature in module.feature
            self.__update(round(memberKME.rx(feature, 1)[0],2))


def write(iteration, kME):
    '''
    writes the eigengene connectivity (kME)
    dictionary to file
    '''
    df = ro.DataFrame(kME)
    df.rownames = (iteration)
    write_data_frame(df, 'eigengene-connectivity.txt', 'Iteration')
