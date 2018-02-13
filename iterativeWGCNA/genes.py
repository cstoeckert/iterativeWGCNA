# pylint: disable=invalid-name
# pylint: disable=unused-import
'''
manage genes
'''
from __future__ import print_function

import logging
from collections import OrderedDict
from collections import Counter

import rpy2.robjects as ro

# from .expression import Expression
from .analysis import calculate_kME
from .eigengenes import Eigengenes
from .r.imports import wgcna, stats, base, rsnippets, grdevices
from .io.utils import xstr
from .r.manager import RManager

class Genes(object):
    '''
    track input genes and their properties, including
    expression profiles, module membership, eigengene
    connectivity
    '''

    def __init__(self, exprData, debug=False):
        '''
        initialize an OrderedDict of genes
        from the row.names of the expression
        data set
        '''
        self.logger = logging.getLogger('iterativeWGCNA.Genes')
        self.profiles = exprData
        self.genes = OrderedDict((geneId, {'module': 'UNCLASSIFIED',
                                           'kME':float('NaN'),
                                           'iteration': None})
                                 for geneId in self.profiles.genes())

        self.size = len(self.genes)
        self.iteration = None
        self.debug = debug


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


    def __update_classified_iteration(self, gene, iteration):
        '''
        set the iteration during which
        a gene was first classified
        '''
        if gene in self.genes:
            self.genes[gene]['iteration'] = iteration
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
                module = self.iteration + '_' + 'M' + str(module)
                self.__update_classified_iteration(gene, self.iteration)
            self.__update_module(gene, module)

        return None


    def copy_membership(self, source):
        '''
        updates membership from another Genes object
        '''
        sourceMembership = source.get_gene_membership()
        for gene, module in sourceMembership.items():
            self.__update_module(gene, module)


    def __extract_iteration_genes(self, targetIteration):
        '''
        get genes classified during specified interation
        '''
        assignedIterations = self.__extract_classified_iteration()
        return [gene for gene, iteration in assignedIterations.items()
                if iteration == targetIteration]


    def __extract_classified_iteration(self):
        '''
        get classified iteration as an ordered dict
        '''
        return OrderedDict((gene, membership['iteration'])
                           for gene, membership in self.genes.items())


    def __extract_modules(self):
        '''
        extract module membership as an ordered dict
        '''
        return OrderedDict((gene, membership['module']) for gene, membership in self.genes.items())


    def get_gene_membership(self, genes=None):
        '''
        public facing method for getting gene membership; returns a
        gene -> membership hash

        if gene list is provided, return only the membership assignment
        for the provided genes
        '''
        if genes is None:
            return self.__extract_modules()
        else:
            return OrderedDict((gene, module) for gene, module in self.__extract_modules().items() if gene in genes)


    def get_gene_kME(self):
        '''
        public facing method for getting all gene kMEs
        '''
        return self.__extract_kME()


    def get_iteration_kME(self, iteration):
        '''
        return kME for all assignments made
        during current iteration
        '''
        iterationGenes = self.__extract_iteration_genes(iteration)
        geneKME = self.__extract_kME()
        return [kME for gene, kME in geneKME.items() if gene in iterationGenes]


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


    def write(self, prefix='', iteration=None):
        '''
        writes the membership and eigengene connectivity
        to files
        filtering for specific iteration if specified
        '''
        summaryGenes = None
        if iteration is None:
            summaryGenes = self.genes
        else:
            iterationGenes = self.__extract_iteration_genes(iteration)
            summaryGenes = OrderedDict((gene, membership) for gene, membership
                                       in self.genes.items()
                                       if gene in iterationGenes
                                       and membership['module'] != 'UNCLASSIFIED')

        with open(prefix + 'membership.txt', 'w') as f:
            print('\t'.join(('Gene', 'Module', 'kME')), file=f)
            for g in summaryGenes:
                print('\t'.join((g, self.genes[g]['module'], xstr(self.genes[g]['kME']))), file=f)
        return None


    def write_iteration_counts(self, prefix=''):
        '''
        print iteration summary
        '''
        with open(prefix + 'summary.txt', 'w') as f:
            print('\t'.join(('N Input Genes', 'N Classified Genes',
                             'N Residual Genes', 'N Detected Modules')), file=f)
            numClassifiedGenes = self.count_classified_genes()
            print('\t'.join((str(self.size),
                             str(numClassifiedGenes),
                             str(self.size - numClassifiedGenes),
                             str(self.count_modules(self.get_classified_genes())))), file=f)


    def plot_kme_histogram(self, iteration, prefix='', vline=0.80):
        '''
        generate kme histogram for genes classified in
        current iteration
        '''
        kmeVector = None
        if 'final' in prefix or 'merge' in prefix:
            classifiedGenes = self.get_classified_genes()
            geneKME = self.__extract_kME()
            kmeVector = [kME for gene, kME in geneKME.items() if gene in classifiedGenes]
        else:
            kmeVector = self.get_iteration_kME(iteration)

        if kmeVector is not None:
            if len(kmeVector) != 0:
                manager = RManager(kmeVector)
                grdevices().pdf(prefix + "kme_histogram.pdf")
                manager.histogram(vline, {'main': 'Gene -> Assigned Module kME for iteration ' + iteration,
                                          'xlab': 'kME', 'ylab':'Gene Count'})
                grdevices().dev_off()


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
        classified = [gene for gene, module in membership.items()
                      if module != 'UNCLASSIFIED']

        return len(classified)



    def get_classified_genes(self, genes=None):
        '''
        gets the list of classified genes
        if a list of genes is provided, only returns
        genes within the specified list
        '''
        membership = self.__extract_modules()
        if genes is not None:
            membership = OrderedDict((gene, module) for gene, module in membership.items()
                                     if gene in genes)
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
                self.__update_classified_iteration(g, None)


    def get_modules(self, genes=None):
        '''
        gets list of unique modules from gene membership assignments
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
                self.__update_classified_iteration(g, None)

            if kME < minKMEtoStay:
                self.__update_module(g, 'UNCLASSIFIED')
                self.__update_kME(g, float('NaN'))
                self.__update_classified_iteration(g, None)


    def merge_close_modules(self, eigengenes, cutHeight):
        '''
        merge close modules based on similarity between
        eigengenes

        return updated eigengene object
        '''

        # repeat until no more merges are possible
        noMergesFound = False
        mergeCount = 0
        modules = self.get_modules()
        classifiedGenes = self.get_classified_genes()
        classifiedGeneProfiles = self.profiles.gene_expression(classifiedGenes)
        while not noMergesFound:
            # compare modules, finding min dissimilarity
            similarity = eigengenes.similarity(None)
            closeModules = rsnippets.findCloseModules(similarity, cutHeight)
            if closeModules != ro.NULL:
                m1 = closeModules.rx2('m1')[0]
                m2 = closeModules.rx2('m2')[0]
                dissimilarity = closeModules.rx2('dissimilarity')[0]
                mergeCount = mergeCount + 1
                self.logger.info("Merging " + m1 + " into " + m2
                                 + " (D = " + str(dissimilarity) + ")")

                memberGenes = self.get_module_members(m1)
                for g in memberGenes:
                    self.__update_module(g, m2)
                    self.__update_classified_iteration(g, 'FINAL_MERGE')
                self.__update_module_kME(m1, eigengenes.get_module_eigengene(m2))

                modules = self.get_modules()
                classifiedGeneMembership = self.get_gene_membership(classifiedGenes)
                if self.debug:
                    self.logger.debug("Getting module assignments for classified genes")
                    self.logger.debug(classifiedGeneMembership)

                eigengenes.recalculate(classifiedGeneProfiles,
                                       classifiedGeneMembership)

            else:
                noMergesFound = True

        self.logger.info("Done merging close modules: " + str(mergeCount) + " modules merged.")
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


    def load_membership(self, fileName=None):
        '''
        loads membership
        '''
        if fileName is None:
            fileName = "final-membership.txt"

        membership = ro.DataFrame.from_csvfile(fileName, sep='\t',
                                               header=True, row_names=1, as_is=True)

        if self.debug:
            self.logger.debug("Loaded membership from file " + fileName + "; see R-log")
            ro.r("print('Loaded membership from file -- head of file:')")
            self.logger.debug(membership.head())

        index = membership.names.index('Module') + 1 # add 1 b/c of python/rpy2/R inconsistency

        if self.debug:
            self.logger.debug("Adjusted index of Module column: " + str(index))

        classifiedCount = 0
        unclassifiedCount = 0
        for g in self.genes:
            gStr = ro.StrVector([str(g)])
            # strange, but necessary so that rpy2 will treat numeric gene ids as strings
            # python str() conversion did not work

            module = membership.rx(gStr[0], index)[0]

            if module == 'UNCLASSIFIED':
                unclassifiedCount = unclassifiedCount + 1
            else:
                classifiedCount = classifiedCount + 1
            self.__update_module(g, module)

        self.logger.info("Loaded " + str(classifiedCount) + " classified genes")
        self.logger.info("Loaded " + str(unclassifiedCount) + " unclassified genes")
