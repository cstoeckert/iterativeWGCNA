'''
a module is a cluster of sequence features
identified by the WGCNA analysis
the module classs links the member sequence
features to their expression profiles,
eigengenes, and eigengene connectivity
'''

from .r.imports import wgcna, stats, base
import iwgcna.membership as m
import iwgcna.eigengene_connectivity as ec

class IterativeWGCNA(object):

    def __init__(self, expr):  
        self.expression = expr
        self.membership = Membership(self.expression)
        self.kme = EigengeneConnectivity(self.expression)
        
    def get_expression(self):
        return self.expression
        
    def calculate_kme(self, eigengene, calculateP):
        '''
        calculates eigengene connectivity (kme)
        between an eigengene and the expression data set
        '''
        if calculateP:
            localKME = wgcna().corAndPvalue(base().t(self.expression), base().t(eigengene))
        else:
            localKME = base().as_data_frame(stats().cor(base().t(self.expression), base().t(eigengene)))

        return localKME

   
