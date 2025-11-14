function [knockoutGrowth, knockoutRatios] = runReactionDeletiononModels(cell_models, binaryRxnsTable, method)

knockoutRatios = binaryRxnsTable;
knockoutGrowth = binaryRxnsTable;
nmodels = length(cell_models);
model_names = binaryRxnsTable.Properties.VariableNames(5:end);
%availableMets = erase(returnAvailableMets(), '[e]');

%% Optimize for growth 
names = cell(nmodels, 1); 
for i = 1:nmodels 
    curModel = cell_models{i}; 
    curModel.csense = repmat('E', length(curModel.mets), 1); 
    
    fprintf('\n====================\n');
    fprintf('Knocking out model %s...',  curModel.id);

    % Only allow media
    %[curModel] = setExchangeBounds(curModel, availableMets, -1000, 1000, true, true); 

    % run FBA with minsum
    %curSol = solveLP(curModel, minFlux);
    %solWT = optimizeCbModel(curModel, 'max', 'one');

    [grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] =singleRxnDeletion(curModel, method);
 
    % Find where the 1s are in the binaryTable
    idx = knockoutRatios.(model_names{i}) == 1;  

    model_name = model_names{i};
    knockoutGrowth.(model_name)(idx) = grRateKO; 
    knockoutRatios.(model_name)(idx) = grRatio;
    
end 