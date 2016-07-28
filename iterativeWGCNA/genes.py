# pylint: disable=invalid-name
# pylint: disable=unused-import
'''
manage genes
'''

import logging
from collections import OrderedDict
from collections import Counter

import rpy2.robjects as ro

# from .expression import Expression
from .analysis import calculate_kME
from .r.imports import wgcna, stats, base, rsnippets
from .io.utils import write_data_frame

class Genes(object):
    '''
    track input genes and their properties, including
    expression profiles, module membership, eigengene
    connectivity
    '''

    def __init__(self, exprData):
        '''
        initialize an OrderedDict of genes
        from the row.names of the expression
        data set
        '''
        self.logger = logging.getLogger('iterativeWGCNA.Genes')
        self.profiles = exprData
        self.genes = OrderedDict((geneId, {'module': 'UNCLASSIFIED', 'kME':float('NaN')}) for geneId in self.profiles.genes())

        self.size = len(self.genes)
        self.iteration = None


    def get_module(self, gene):
        '''
        returns the assigned module for a gene
        '''
        return self.genes[gene]['module']


    def __is_classified(self, gene):
        '''
        returns true if the feature is classified
        '''
        return self.get_module(gene) != 'UNCLASSIFIED'


    def __update_module(self, gene, module):
        '''
        update gene module
        do not add new genes
        '''
        if gene in self.genes:
            self.genes[gene]['module'] = module
            return True
        else:
            return False


    def update_membership(self, genes, blocks):
        '''
        fetches new module membership from WGCNA
        blocks and updates relevant genes
        '''
        modules = rsnippets.extractModules(blocks, genes)
        # if the feature is in the subset
        # update, otherwise leave as is
        for gene in genes:
            # .rx returns a FloatVector which introduces
            # a .0 to the numeric labels when converted to string
            # which needs to be removed
            # note: R array starts at index 1, python at 0
            module = str(modules.rx(gene, 1)[0]).replace('.0', '')
            if module in ('0', 'grey'):
                module = 'UNCLASSIFIED'
            else:
                module = self.iteration + '-' + module
            self.__update_module(gene, module)

        return None


    def copy_membership(self, source):
        '''
        updates membership from another Genes object
        '''
        sourceMembership = source.get_gene_membership()
        for gene, module in sourceMembership.items():
            self.__update_module(gene, module)
        
    
    def __extract_modules(self):
        '''
        extract module membership as an ordered dict
        '''
        return OrderedDict((gene, membership['module']) for gene, membership in self.genes.items())


    def get_gene_membership(self):
        '''
        public facing method for getting gene membership
        '''
        return self.__extract_modules()

    def __extract_kME(self):
        '''
        extract eigengene connectivity (kME)
        as an ordered dict
        '''
        return OrderedDict((gene, membership['kME']) for gene, membership in self.genes.items())


    def __write_modules(self, isPruned):
        '''
        writes the gene memebership to a file
        '''
        df = ro.DataFrame(self.__extract_modules())
        df.rownames = (self.iteration)
        fileName = 'membership.txt' if isPruned else 'pre-pruning-membership.txt'
        write_data_frame(df, fileName, 'Iteration')
        return None


    def __update_kME(self, gene, kME):
        '''
        update gene eigengene connectivity (kME)
        do not add new genes
        '''
        if gene in self.genes:
            self.genes[gene]['kME'] = kME
            return True
        else:
            return False


    def __update_module_kME(self, module):
        '''
        update member gene eigengene connectivity (kME)
        for specified module (module is an object
        containing module name and eigengene)
        '''

        members = self.get_module_members(module.name)
        memberKME = calculate_kME(self.profiles.get_expression(members),
                                  module.eigengene, False)

        self.logger.debug("updating module kme")
        self.logger.debug("Module Name:" +  module.name)
        self.logger.debug("Module members: " + members)
        self.logger.debug(memberKME)

        for gene in memberKME.rownames:
            self.__update_kME(gene, round(memberKME.rx(gene, 1)[0], 2))

    def update_kME(self):
        modules = self.get_modules()
        for m in modules:
            self.__update_module_kME(m)

    def __write_kME(self):
        '''
        writes eigengene connectivity (kME)
        to a file
        '''
        df = ro.DataFrame(self.__extract_kME())
        df.rownames = (self.iteration)
        fileName = 'eigengene-connectivity.txt'
        write_data_frame(df, fileName, 'Iteration')
        return None


    def write(self, isPruned):
        '''
        writes the membership and eigengene connectivity
        to files
        '''
        self.__write_modules(isPruned)
        self.__write_kME()
        return None


    def count_module_members(self, genes=None):
        '''
        counts the number of genes per module
        and returns a dict of module -> gene count
        if a list of genes is provided, only counts within
        the specified gene list
        '''
        membership = self.__extract_modules()
        if genes is not None:
            membership = {gene:module for gene, module in membership.items() if gene in genes}
        return Counter(membership.values())


    def count_classified_genes(self, genes=None):
        '''
        counts and return the number of classified genes
        if a list of genes is provided, only counts within
        the specified gene list
        '''

        membership = self.__extract_modules()
        if genes is not None:
            membership = {gene:module for gene, module in membership.items() if gene in genes}
        classified = {gene:module for gene, module in membership.items()
                      if module != 'UNCLASSIFIED'}

        return len(classified)


    def get_classified_gene_profiles(self):
        '''
        get classified gene profiles
        '''
        membership = self.__extract_modules()
        classifiedGenes = [gene for gene, module in membership.items() if module != 'UNCLASSIFIED']
        return self.profiles.get_expression(classifiedGenes)


    def get_unclassified_gene_profiles(self):
        '''
        get unclassified gene profiles
        '''
        membership = self.__extract_modules()
        unclassifiedGenes = [gene for gene, module in membership.items() if module == 'UNCLASSIFIED']
        return self.profiles.get_expression(unclassifiedGenes)


    def get_gene_expression(self, genes):
        '''
        get expression profiles for specified
        list of genes
        '''
        return self.profiles.get_expression(genes)

    
    def count_modules(self, genes=None):
        '''
        counts the number of modules (excluding unclassified)
        if a list of genes is provided, only counts within
        the specified gene list
        '''
        moduleCount = self.count_module_members(genes)
        return len(moduleCount) - 1 if 'UNCLASSIFIED' in moduleCount else len(moduleCount)


    def remove_small_modules(self, minModuleSize):
        '''
        checks membership counts and removes
        any modules that are too small
        by updating gene membership to UNCLASSIFIED and
        setting eigengene connectivity (kME) to NaN
        '''
        memberCount = self.count_module_members()

        for gene, module in self.genes.items():
            if memberCount[module] < minModuleSize:
                self.__update_module(gene, 'UNCLASSIFIED')
                self.__update_kME(gene, float('NaN'))


    def get_modules(self):
        '''
        gets list of modules from gene membership assignments
        '''
        # get unique members by converting values to a set
        membership = set(self.__extract_modules().values())
        membership.discard('UNCLASSIFIED')
        return list(membership)


    def get_module_members(self, targetModule):
        '''
        get list of module member genes
        '''
        membership = self.__extract_modules()
        return [gene for gene, module in membership.items() if module == targetModule]
