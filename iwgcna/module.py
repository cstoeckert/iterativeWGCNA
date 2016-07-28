'''
a module is a cluster of sequence features
identified by the WGCNA analysis
the module classs links the member sequence
features to their expression profiles,
eigengenes, and eigengene connectivity
'''

from .r.imports import wgcna, stats, base

class Module(object):

    def __init__(self, name, features, expr):
        self.name = name
        self.features = features
        self.expression = expr
        self.eigengene = None
        self.kme = None

    def calculate_kme(self, calculateP):
        '''
        calculates eigengene connectivity kme
        between an eigengene and expression data set
        '''
        if calculateP:
            self.kme = wgcna().corAndPvalue(base().t(self.expr), base().t(self.eigengene))
        else:
            self.kme = base().as_data_frame(stats().cor(base().t(self.expr), base().t(self.eigengene)))

        return None

    def get_kme(self):
        '''
        return module kme (eigengene connectivity)
        between every member feature and the module
        eigengene
        '''
        return self.kme
