'''
manage module membership lists
'''

import logging
from collections import OrderedDict
import rpy2.robjects as ro

from .r.imports import r_utils
from .io.utils import write_data_frame
from . import kme


class Membership(object):
    ''' track module membership for the sequence features in the expression dataset'''

    def __init__(self, features):
        '''
        initialize an OrderedDict with one entry per feature
        all features are initially unclassified
        '''
        self.classification = OrderedDict((f, 'UNCLASSIFIED') for f in features)
        self.size = len(features)
        self.iteration = None
        return None

    def __update(self, feature, module):
        '''
        update classification of 'feature' to 'module'
        do not add new features
        '''
        if feature in self.classification:
            self.classification[feature] = module
            return True
        else:
            return False

    def set_current_iteration(self, iteration):
        ''' sets current iteration '''
        self.iteration = iteration
        return None

    def update(self, features, blocks):
        '''
        compares new module membership assignments to
        prexisting ones; updates classification
        '''
        modules = r_utils.extractModules(blocks, features)
        # if the feature is in the subset
        # update, otherwise leave as is
        for f in features:
            # .rx returns a FloatVector which introduces
            # a .0 to the numeric labels when converted to string
            # which needs to be removed
            # note: R array starts at index 1, python at 0
            module = str(modules.rx(g, 1)[0]).replace('.0', '')
            if module in ('0', 'grey'):
                module = 'UNCLASSIFIED'
            else:
                module = self.iteration + '-' + module
            self.__update(f, module)

        return None

    def write(self, isPruned):
        '''
        writes the membership dictionary to file
        :param iteration    iWGCNA iteratoin
        :param membership   gene->module mapping
        :param initialClassificaton    boolean flag indicating whether pruning has been done
        '''
        df = ro.DataFrame(self.classification)
        df.rownames = (self.iteration)
        fileName = 'membership.txt' if isPruned else 'pre-pruning-membership.txt'
        write_data_frame(df, fileName, 'Iteration')

        return None

    def getClassification(self, feature):
        '''
        returns the assigned module for a feature
        '''
        return self.classification[feature]

    def count_module_members(self, features=None):
        '''
        counts the number of features per module
        and returns a dict of module -> feature count
        if a list of features is provided, only counts within
        the specified feature list
        '''
        if features is None:
            features = list(self.classification.keys())

        count = {}
        for feature in features:
            module = self.getClassification(feature)
        if module in count:
            count[module] = count[module] + 1
        else:
            count[module] = 1

        return count

    def isClassified(self, feature):
        '''
        returns true if the feature is classified
        '''

        return self.getClassification(feature) != 'UNCLASSIFIED'

    def count_classified_features(self, features=None):
        '''
        counts and return the number of classified features
        if a list of features is provided, only counts within
        the specified feature list
        '''
        count = 0
        if features is None:
            features = list(self.classification.keys())
        for feature in features:
            if self.isClassified(feature):
                count = count + 1

        return count

    def count_modules(self, features=None):
        '''
        counts the number of modules (excluding unclassified)
        if a list of features is provided, only counts within
        the specified feature list
        '''
        moduleCount = self.count_module_members(features)
        return len(moduleCount) - 1 if 'UNCLASSIFIED' in moduleCount else len(moduleCount)

    def remove_small_modules(self, kME, minModuleSize):
        '''
        checks membership counts and removes
        any modules that are too small
        by updating feature membership to UNCLASSIFIED and setting kME to NaN
        returns update kME and membership
        '''
        memberCount = self.count_module_members()

        for feature, module in self.classification.items():
            if memberCount[module] < minModuleSize:
                self.__update(feature, 'UNCLASSIFIED')
                kME[feature] = float('NaN')
                # TODO: kME.update(feature, float('NaN'))
                # move update of the kME out of here?

        return kME

    def get_modules(self):
        '''
        gets list of modules in  membership assignments
        '''
        modules = {}
        for module in self.classification.values():
            if module != "UNCLASSIFIED":
                modules[module] = 1

        return list(modules.keys())


    def get_members(self, target):
        '''
        get list of features in specified module
        '''
        members = []
        for feature, module in self.classification.items():
            if module == target:
                members.append(feature)

        return members


    # TODO -- move this out of membership and into iwgcna

    def best_fit(membership, eigengenes, data, kME, params):
        '''
        Evaluate eigengene connectivity (kME)
        for each feature against the eigengenes for each
        of the final modules found by iWGCNA.
        If kME(module) > kME(assigned_module)
        and the p-value <= the reassignThreshold (of WGCNA
        parameters) then reassign the module
        membership of the feature.
        '''
        reassignmentCount = 0

        for module in eigengenes.rownames:
            logging.info("Evaluating best fits to " + module)
            # calculate kME of all features to the module eigengene
            moduleEigengene = eigengenes.rx(module, True)
            moduleKME = kme.calculate(data, moduleEigengene, True)
            
            # evaluate each feature
            for feature in membership:
                currentModule = membership[feature]

                # I believe this check is no longer necessary b/c we
                # should have already filtered out duplicate eigengenes
                # from earlier iterations
                # if module == currentModule or \
                #     eigen_equal(moduleEigengene, eigengenes.rx(currentModule, True)):
                if module == currentModule:
                    continue
            
                currentKME = kME[feature]

                newKME = round(moduleKME.rx2('cor').rx(feature, 1)[0], 2)
                pvalue = moduleKME.rx2('p').rx(feature, 1)[0]

                if (currentModule == "UNCLASSIFIED" \
                    and newKME >= params['minKMEtoStay']) \
                    or (newKME > currentKME \
                    and pvalue < params['reassignThreshold']):
        
                    membership[feature] = module
                    kME[feature] = newKME
                    reassignmentCount = reassignmentCount + 1

        return reassignmentCount, kME

