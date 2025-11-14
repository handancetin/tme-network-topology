function [maxGrowths, solutionsTable] = runMinFBAonModels(cell_models, binaryRxnsTable, minFlux)

solutionsTable = binaryRxnsTable;
nmodels = length(cell_models);
model_names = binaryRxnsTable.Properties.VariableNames(5:end);
%availableMets = erase(returnAvailableMets(), '[e]');

%% Optimize for growth
maxGrowths = zeros(nmodels, 1); 
names = cell(nmodels, 1); 
for i = 1:nmodels 
    curModel = cell_models{i}; 
    
    fprintf('\n====================\n');
    fprintf('Simulating model %s...',  curModel.id);

    % Only allow media
    %[curModel] = setExchangeBounds(curModel, availableMets, -1000, 1000, true, true); 

    % run FBA with minsum
    curSol = solveLP(curModel, minFlux);
 
    % Find where the 1s are in the binaryTable
    idx = solutionsTable.(model_names{i}) == 1; 

    model_name = model_names{i};
    solutionsTable.(model_name)(idx) = curSol.x;

    maxGrowths(i) = curSol.f;

end
  