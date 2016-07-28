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

