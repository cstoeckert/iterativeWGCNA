# pylint: disable=invalid-name
# pylint: disable=unused-import
'''
stores a single module; for summary purposes
'''

import logging

import rpy2.robjects as ro

# from .expression import Expression

# from .eigengenes import Eigengenes
# from .r.imports import wgcna, stats, base, rsnippets
# from .io.utils import write_data_frame

from .analysis import calculate_kME

class Module(object):
    '''
    a module, incl: size, members, expression profiles, name, color
    used to facilitate per module summaries
    '''

    def __init__(self, name, genes, eigengene):
        '''
        initialize module object
        '''
        self.logger = logging.getLogger('iterativeWGCNA.Module')
        self.name = name
        self.members = genes.get_module_members(self.name)
        self.expression = None
        self.size = len(self.members)
        self.eigengene = eigengene
        self.kME = genes.get_module_kME(self.name)


    def plot_eigengene(self):
        ''' 
        plot module eigengene
        '''
        # barplot(x, xlim=c(-1,1), col=as.vector(ifelse(data[2,]>0,"red","blue")), horiz=T, names.arg=names(data), main = "Eigengene", ylab="Sample", xlab="Standardized Expression", las=1, cex.names =0.5)
