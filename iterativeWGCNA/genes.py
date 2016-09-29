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
from .eigengenes import Eigengenes
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
        modules = rsnippets.extractModules(blocks, ro.StrVector(genes))
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


    def get_gene_kME(self):
        '''
        public facing method for getting all gene kMEs
        '''
        return self.__extract_kME()


    def get_module_kME(self, targetModule):
        '''
        get all kME values in a module
        '''
        membership = self.get_module_members(targetModule)
        memberKME = self.__extract_kME()
        return [kME for gene, kME in memberKME.items() if gene in membership]


    def get_kME(self, gene):
        '''
        returns the assigned kME for a gene
        '''
        return self.genes[gene]['kME']


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


    def __update_module_kME(self, module, eigengene, genes=None):
        '''
        update member gene eigengene connectivity (kME)
        for specified module and eigengene
        '''
        members = self.get_module_members(module)
        memberKME = calculate_kME(self.profiles.gene_expression(members),
                                  eigengene, False)

        for gene in memberKME.rownames:
            if genes is not None:
                if gene in genes:
                    self.__update_kME(gene, round(memberKME.rx(gene, 1)[0], 2))


    def update_kME(self, eigengenes, genes=None):
        '''
        update module kME given its eigengene
        '''
        modules = self.get_modules()
        for m in modules:
            moduleEigengene = eigengenes.get_module_eigengene(m)
            self.__update_module_kME(m, moduleEigengene, genes)


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
        if isPruned:
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



    def get_classified_genes(self, genes=None):
        '''
        gets the list of classifed genes
        if a list of genes is provided, only returns
        genes within the specified list
        '''
        membership = self.__extract_modules()
        if genes is not None:
            membership = {gene:module for gene, module in membership.items()
                          if gene in genes}
        classifiedGenes = [gene for gene, module in membership.items()
                           if module != 'UNCLASSIFIED']

        return classifiedGenes



    def get_unclassified_genes(self):
        '''
        get unclassified genes
        '''
        membership = self.__extract_modules()
        unclassifiedGenes = [gene for gene, module in membership.items()
                             if module == 'UNCLASSIFIED']
        return unclassifiedGenes



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

        for g in self.genes:
            geneModule = self.get_module(g)
            if memberCount[geneModule] < minModuleSize:
                self.__update_module(g, 'UNCLASSIFIED')
                self.__update_kME(g, float('NaN'))


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


    def get_genes(self):
        '''
        return list of all genes
        '''
        return [gene for gene in self.genes]


    def evaluate_fit(self, minKMEtoStay, genes=None):
        '''
        evaluate fit of each gene to its assigned
        module, unclassifying if the fit is below the
        minimum KME to stay
        if a gene list is provided, only evaluates the
        specified genes
        '''

        if genes is None:
            genes = self.profiles.genes()

        for g in genes:
            module = self.get_module(g)
            kME = self.get_kME(g)

            if module == 'UNCLASSIFIED':
                self.__update_kME(g, float('NaN'))

            if kME < minKMEtoStay:
                self.__update_module(g, 'UNCLASSIFIED')
                self.__update_kME(g, float('NaN'))


    def merge_close_modules(self, eigengenes, cutHeight):
        '''
        merge close modules based on similarity between
        eigengenes

        return updated eigengene object
        '''
        modules = self.get_modules()
        revisedModules = {}

        for m1 in modules:
            similarity = eigengenes.similarity(m1)
            for m2 in modules:
                if m1 == m2:
                    continue
                dissimilarity = round(1.0 - similarity.rx(m2, 1)[0], 2)
                if dissimilarity <= cutHeight:
                    revisedModules[m1] = m2

                    # remove m2 so we don't end up mapping m1 to m2 as well as m2 to m1
                    modules.remove(m2)

                    self.logger.info("Merging " + m1 + " and "
                                     + m2 + " (D = " + str(dissimilarity) + ")")

        if len(revisedModules) != 0:
            # update gene membership and kME calculations
            for m in revisedModules:
                memberGenes = self.get_module_members(m)
                newModule = revisedModules[m]
                for g in memberGenes:
                    self.__update_module(g, newModule)
                self.__update_module_kME(m, eigengenes.get_module_eigengene(newModule))

                modules = self.get_modules() # update list of modules
                eigengenes.update_to_subset(modules) # update eigengenes to reflect new list

        self.logger.info("Done merging close modules: " + str(len(revisedModules)) + " modules merged.")
        self.logger.info("Retained " + str(len(modules)) + " modules after merge.")
        return eigengenes


    def reassign_to_best_fit(self, eigengenes, reassignThreshold, minKMEtoStay):
        '''
        Evaluate eigengene connectivity (kME)
        for each feature against the eigengenes for each
        of the final modules found by iWGCNA.
        If kME(module) > kME(assigned_module)
        and the p-value <= the reassignThreshold (of WGCNA
        parameters) then reassign the module
        membership of the feature.

        returns a count of the number of reassigned genes
        '''
        count = 0
        modules = self.get_modules()
        for m in modules:
            # calculate kME of all genes to the module eigengene
            moduleEigengene = eigengenes.get_module_eigengene(m)
            moduleKME = calculate_kME(self.profiles.expression(), moduleEigengene, True)

            # for each gene not assigned to the current module, test fit
            for g in self.genes:
                currentModule = self.get_module(g)
                if currentModule != m:
                    kME = self.get_kME(g)
                    newKME = round(moduleKME.rx2('cor').rx(g, 1)[0], 2)
                    pvalue = moduleKME.rx2('p').rx(g, 1)[0]

                    if (currentModule == "UNCLASSIFIED" \
                        and newKME >= minKMEtoStay) \
                        or (newKME > kME \
                        and pvalue < reassignThreshold):

                        self.__update_module(g, m)
                        self.__update_kME(g, newKME)
                        count = count + 1

        return count

