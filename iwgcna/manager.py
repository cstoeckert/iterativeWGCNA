'''
functions for running iterative WGCNA
'''
from __future__ import print_function

import logging

from .r.imports import base, wgcna, r_utils
import iwgcna.eigengenes as eigengenes
import iwgcna.membership as membership
import iwgcna.kme as kme

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
        membershipMap, kmeMap = \
          evaluate_fit(kmeMap, membershipMap, data.rownames,
                       wgcnaParameters['minKMEtoStay'])
        membershipMap, kmeMap = \
          membership.remove_small_modules(membershipMap, kmeMap,
                                          wgcnaParameters['minModuleSize'])
        membership.write(iteration, membershipMap, True)

    return membershipMap, kmeMap


def merge_close_modules(eigengeneMatrix, membershipMap, cutHeight):
    '''
    merge close modules based on similarity between
    eigengenes
    '''
    revisedModules = {}
    modules = membership.get_modules(membershipMap)
    for module1 in modules:
        similarity = eigengenes.similarity(eigengeneMatrix, module1)
        for module2 in modules:
            if module1 == module2:
                continue
            dissimilarity = round(1.0 - similarity.rx(module2, 1)[0], 2)
            if dissimilarity <= cutHeight:
                logging.info("Merging " + module1 + " and "
                             + module2 + " (D = " + str(dissimilarity) + ")")

                revisedModules[module1] = module2
                # remove m2 so it is not considered later
                # on down the line and we don't end up
                # mapping m1 to m2 as well as m2 to m1
                modules.remove(module2)

    if len(revisedModules) == 0:
        logging.info("Merge Close Modules: 0 modules merged.")
        return eigengeneMatrix, membershipMap, None

    for gene in membershipMap:
        module = membershipMap[gene]
        if module in revisedModules:
            membershipMap[gene] = revisedModules[module]

    newEigengenes = eigengenes.extract_modules(eigengeneMatrix,
                                               membership.get_modules(membershipMap))
    numMergedModules = eigengeneMatrix.nrow - newEigengenes.nrow
    logging.info("Merge Close Module: " + str(numMergedModules) + " modules merged.")

    return newEigengenes, membershipMap, numMergedModules
