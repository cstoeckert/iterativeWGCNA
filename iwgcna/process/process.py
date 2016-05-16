'''
functions for running iterative WGCNA
'''

from __future__ import print_function
import os
from ..r.imports import base, wgcna, r_utils
import iwgcna.process.eigengenes as eigengenes
import iwgcna.process.membership as membership
import iwgcna.process.kme as kme

def write_gene_counts(iteration, initial, fit):
    '''
    writes the number of kept and dropped genes at the end of an iteration
    '''
    fileName = 'gene-counts.txt'
    try:
        os.stat(fileName)
    except OSError:
        header = ('Iteration', 'Initial', 'Fit', 'Residual')
        with open(fileName, 'a') as f:
            print('\t'.join(header), file=f)
    finally:
        with open(fileName, 'a') as f:
            print('\t'.join((iteration, str(initial),
                             str(fit), str(initial - fit))), file=f)


def set_iteration_label(runId, passId):
    '''
    generates the unique label for the iteration
    '''
    label = 'p' + str(passId) + '_iter_' + str(runId)
    return label


def evaluate_fit(kmeMap, membershipMap, genes, minKMEtoStay):
    '''
    evaluate eigengene similarity (KME) as
    a measure of goodness of fit
    label poorly fitting genes as unclassified
    '''
    for g in genes:
        module = membershipMap[g]

        if module == 'UNCLASSIFIED':
            kmeMap[g] = float('NaN') 
            continue

        if kmeMap[g] < minKMEtoStay:
            membershipMap[g] = 'UNCLASSIFIED'
            kmeMap[g] = float('NaN')

    return membershipMap, kmeMap


def run_wgcna(data, iteration, params):
    '''
    run wgcna
    :param data  expression data
    :param iteration    unique id for the iteration
    :param params    dict storing WGCNA parameter name:value pairs
    '''

    params['datExpr'] = base().t(data) # have to transpose before passing to WGCNA
    params['saveTOMFileBase'] = iteration + '-TOM'

    return wgcna().blockwiseModules(**params)


def run_iteration(iteration, data, membershipMap, kmeMap,
                  wgcnaParameters, saveBlocks):
    '''
    run iteration of blockwise WGCNA
    return membership and eigengenes
    '''

    # run blockwise WGCNA
    blocks = run_wgcna(data, iteration, wgcnaParameters)
    if saveBlocks:
        r_utils.saveObject(blocks, 'blocks', iteration + '-blocks.RData')

    # extract eigengenes from blocks
    # if there are eigengenes, then evaluate fitness
    eigengeneMatrix = r_utils.extractEigengenes(iteration, blocks, data.colnames)
    if eigengeneMatrix.nrow != 0:
        eigengenes.write(eigengeneMatrix, False)

        # extract module membership from blocks and update
        membershipMap = membership.update(iteration, data.rownames, blocks, membershipMap)
        membership.write(iteration, membershipMap, False)

        # calculate eigengene connectivity (kME) to
        # assigned module
        kmeMap = kme.update(kmeMap, data, membershipMap, eigengeneMatrix)
        kme.write(iteration, kmeMap)

        # evaluate fit & update membership again, removing small modules
        membershipMap, kmeMap = evaluate_fit(kmeMap, membershipMap, data.rownames,
                                             wgcnaParameters['minKMEtoStay'])
        membershipMap, kmeMap = membership.remove_small_modules(membershipMap, kmeMap,
                                                                wgcnaParameters['minModuleSize'])
        membership.write(iteration, membershipMap, True)

    return membershipMap, kmeMap
