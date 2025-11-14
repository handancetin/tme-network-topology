function [knockoutGrowth, knockoutRatios] = runGeneDeletionOnModels(cell_models, binaryGenesTable, method)

knockoutRatios = binaryGenesTable;
knockoutGrowth = binaryGenesTable;
nmodels = length(cell_models);
model_names = binaryGenesTable.Properties.VariableNames(4:end);
%availableMets = erase(returnAvailableMets(), '[e]');

%% Optimize for growth 
names = cell(nmodels, 1); 
for i = 1:nmodels 
    curModel = cell_models{i}; 
    fprintf('\n====================\n');
    fprintf('Knocking out model %s...',  curModel.id);

    curModel.csense = repmat('E', length(curModel.mets), 1);
    curModel.rules=grrulesToRules(curModel);
    

    % Only allow media
    %[curModel] = setExchangeBounds(curModel, availableMets, -1000, 1000, true, true); 

    % run FBA with minsum
    %curSol = solveLP(curModel, minFlux);
    %solWT = optimizeCbModel(curModel, 'max', 'one');

    [grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] =singleGeneDeletion(curModel, method);
 
    % Find where the 1s are in the binaryTable
    idx = knockoutRatios.(model_names{i}) == 1;  

    model_name = model_names{i};
    knockoutGrowth.(model_name)(idx) = grRateKO; 
    knockoutRatios.(model_name)(idx) = grRatio;
    
end 
end 

function rules=grrulesToRules(model)
    %This function just takes grRules, changes all gene names to
    %'x(geneNumber)' and also changes 'or' and 'and' relations to corresponding
    %symbols
    replacingGenes=cell([size(model.genes,1) 1]);
    for i=1:numel(replacingGenes)
        replacingGenes{i}=strcat('x(',num2str(i),')');
    end
    rules = strcat({' '},model.grRules,{' '});
    for i=1:length(model.genes)
        rules=regexprep(rules,[' ' model.genes{i} ' '],[' ' replacingGenes{i} ' ']);
        rules=regexprep(rules,['(' model.genes{i} ' '],['(' replacingGenes{i} ' ']);
        rules=regexprep(rules,[' ' model.genes{i} ')'],[' ' replacingGenes{i} ')']);
    end
    rules=regexprep(rules,' and ',' & ');
    rules=regexprep(rules,' or ',' | ');
    rules=strtrim(rules);
end