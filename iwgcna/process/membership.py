'''
manage module membership lists
'''

import logging
from collections import OrderedDict
import rpy2.robjects as ro
from ..r.imports import r_utils
from ..utils.io import write_data_frame, warning
# from .eigengenes import equal as eigen_equal
from . import kme

def initialize(data):
    '''
    initialize membership dictionary
    gene list comes from input data row names (DATA.rownames)
    all genes are initially unclassified
    an ordered dictionary is used to keep values in same order
    as input data
    '''

    membership = OrderedDict((gene, 'UNCLASSIFIED') for gene in data.rownames)
    return membership


def update(iteration, genes, blocks, membership):
    '''
    compares new module membership assignments to
    prexisting ones; updates membership list
    '''

    modules = r_utils.extractModules(blocks, genes)

    # if the gene is in the subset
    # update, otherwise leave as is
    for g in genes:
        # .rx returns a FloatVector which introduces
        # a .0 to the numeric labels when converted to string
        # which needs to be removed
        # note: R array starts at index 1, python at 0
        module = str(modules.rx(g, 1)[0]).replace('.0', '')
        if module in ('0', 'grey'):
            module = 'UNCLASSIFIED'
        else:
            module = iteration + '-' + module

        membership[g] = module

    return membership


def write(iteration, membership, isPruned):
    '''
    writes the membership dictionary to file
    :param iteration    iWGCNA iteratoin
    :param membership   gene->module mapping
    :param initialClassificaton    boolean flag indicating whether pruning has been done
    '''
    df = ro.DataFrame(membership)
    df.rownames = (iteration)
    fileName = 'membership.txt' if isPruned else 'pre-pruning-membership.txt'

    write_data_frame(df, fileName, 'Iteration')


def count_module_members(membership, genes=None):
    '''
    counts the number of genes per module
    and returns a dict of module -> gene count
    if a list of genes is provided, only counts within
    the specified gene list
    '''
    if genes is None:
        genes = list(membership.keys())

    count = {}
    for gene in genes:
        module = membership[gene]
        if module in count:
            count[module] = count[module] + 1
        else:
            count[module] = 1

    return count

def count_classified_genes(membership, genes=None):
    '''
    counts and return the number of classified genes
    if a list of genes is provided, only counts within
    the specified gene list
    '''
    count = 0
    if genes is None:
        genes = list(membership.keys())
    for gene in genes:
        if membership[gene] != 'UNCLASSIFIED':
            count = count + 1

    return count


def count_modules(membership, genes=None):
    '''
    counts the number of modules (excluding unclassified)
    if a list of genes is provided, only counts within
    the specified gene list
    '''
    moduleCount = count_module_members(membership, genes)
    return len(moduleCount) - 1 if 'UNCLASSIFIED' in moduleCount else len(moduleCount)


def remove_small_modules(membership, kME, minModuleSize):
    '''
    checks membership counts and removes
    any modules that are too small
    by updating gene membership to UNCLASSIFIED and setting kME to NaN
    returns update kME and membership
    '''
    memberCount = count_module_members(membership)

    for gene in membership:
        module = membership[gene]
        if memberCount[module] < minModuleSize:
            membership[gene] = 'UNCLASSIFIED'
            kME[gene] = float('NaN')
                    
    return membership, kME


def best_fit(membership, eigengenes, data, kME, params):
    '''
    Evaluate eigengene connectivity (kME)
    for each gene against the eigengenes for each
    of the final modules found by iWGCNA.
    If kME(module) > kME(assigned_module)
    and the p-value <= the reassignThreshold (of WGCNA
    parameters) then reassign the module
    membership of the gene.
    '''
    reassignmentCount = 0

    for module in eigengenes.rownames:
        logging.debug("Evaluating best fits to " + module)
        # calculate kME of all genes to the module eigengene
        moduleEigengene = eigengenes.rx(module, True)
        moduleKME = kme.calculate(data, moduleEigengene, True)

        # evaluate each gene
        for gene in membership:
            currentModule = membership[gene]

            # I believe this check is no longer necessary b/c we
            # should have already filtered out duplicate eigengenes
            # from earlier iterations
            # if module == currentModule or \
            #     eigen_equal(moduleEigengene, eigengenes.rx(currentModule, True)):
            if module == currentModule:
                continue

            currentKME = kME[gene]

            newKME = round(moduleKME.rx2('cor').rx(gene, 1)[0], 2)
            pvalue = moduleKME.rx2('p').rx(gene, 1)[0]

            if (currentModule == "UNCLASSIFIED" \
                and newKME >= params['minKMEtoStay']) \
                or (newKME > currentKME \
                and pvalue < params['reassignThreshold']):
      
                membership[gene] = module
                kME[gene] = newKME
                reassignmentCount = reassignmentCount + 1

    return reassignmentCount, membership, kME


def get_modules(membership):
    '''
    gets list of modules in  membership assignments
    '''
    modules = {}
    for gene in membership:
        module = membership[gene]
        if module != "UNCLASSIFIED":
            modules[module] = 1

    return list(modules.keys())


def get_members(module, membership):
    '''
    get list of genes in specified module
    '''
    members = []
    for gene in membership:
        if membership[gene] == module:
            members.append(gene)

    return members
