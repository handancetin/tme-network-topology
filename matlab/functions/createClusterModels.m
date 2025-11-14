function [model, details] = createClusterModels(geneNames, clusterName, mainType, tissue, geneLevels, threshold, prepData)

    fprintf('\n====================\n');
    fprintf('Creating model for: %s...',  clusterName{1});

    % This structure is needed as an input
    transcrData.genes = geneNames;  
    transcrData.tissues = clusterName;
    transcrData.levels = geneLevels;  
    transcrData.threshold =  quantile(geneLevels, threshold);

    % Create the model using ftINIT
     [model, metProduction, addedRxnsForTasks, deletedRxnsInINIT, fullMipRes] = ftINIT(prepData, transcrData.tissues{1}, [], [], transcrData, {}, getHumanGEMINITSteps('1+0'), true, true);
     model.id = transcrData.tissues{1};
     model.maintype = mainType{1};
     model.tissue = tissue{1};
     
     details.metProduction = metProduction;
     details.addedRxnsForTasks = addedRxnsForTasks;
     details.deletedRxnsInINIT = deletedRxnsInINIT;
     details.fullMipRes = fullMipRes;

end