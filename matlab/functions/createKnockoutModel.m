function newmodel = createKnockoutModel(model, target)
    fprintf('Creating a knockout model from %s...\n', model.id)

    newmodel = model;

    % Find reactions that 
    rxnIDs = find(contains(model.eccodes,target));
    fprintf('Total of %d reactions (catalyzed by %s) identified.\n', length(rxnIDs), target)
    
    newmodel.lb(rxnIDs) = 0;
    newmodel.ub(rxnIDs) = 0;
    
    orgsol = solveLP(model);
    sol = solveLP(newmodel);
    if sol.stat
        fprintf('Biomass original %0.2f, biomass new %0.2f, ratio %0.2f.\n', orgsol.f, sol.f, sol.f/orgsol.f)
        newmodel.id = ["Eng" + model.id];
        %sampleUsingRavenSampler({newmodel}, 2000); % done
        %sampleUsingRHMC({newmodel}, 200, 2000) % done
        %saveFilteredRHMCSamples(newmodel);
         createMetabolicWeightedGraphsFromSamples(newmodel);
    else
         fprintf('Knocked out model cannot produce solution. Cancelling.')
         return
    end
    
end