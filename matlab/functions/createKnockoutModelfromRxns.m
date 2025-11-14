function newmodel = createKnockoutModelfromRxns(model, rxns)
    fprintf('Creating a knockout model from %s...\n', model.id)

    newmodel = model;
    newmodel.csense = repmat('E',length(newmodel.mets), 1);
    % Find reactions that 
    rxnIDs = findRxnIDs(model, rxns);
    fprintf('Total of %d reactions will be knocked out.\n', length(rxnIDs))
    
    [grRatio, grRateKO, grRateWT, ~, ~, ~] = singleRxnDeletion(newmodel, 'FBA', rxns, 0);
    
    rxns(grRatio == 0) = [];

    newmodel = removeRxns(newmodel, rxns);
    orgsol = solveLP(model);
    sol = solveLP(newmodel);
    if sol.stat
        fprintf('Biomass original %0.2f, biomass new %0.2f, ratio %0.2f.\n', orgsol.f, sol.f, sol.f/orgsol.f)
        
        newmodel.id = ["Eng" + model.id];
        %sampleUsingRavenSampler({newmodel}, 2000); % done
        %sampleUsingRHMC({newmodel}, 200, 2000) % done
        %saveFilteredRHMCSamples(newmodel);
        %createMetabolicWeightedGraphsFromSamples(newmodel);
    else
         fprintf('Knocked out model cannot produce solution. Cancelling.')
         return
    end
    
end