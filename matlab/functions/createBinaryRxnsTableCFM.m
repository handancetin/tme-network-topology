
function [binaryRxnsMatrix, binaryRxnsTable] = createBinaryRxnsTableCFM(ref_model, subset_models)

    num_models = length(subset_models);
    
    % Initialize the existence matrix
    binaryRxnsMatrix = zeros(length(ref_model.rxns), num_models);
    
    % Extract model IDs and ensure uniqueness
    model_ids = cell(1, num_models);
    for i = 1:num_models
        model_ids{i} = subset_models{i}.id; % Extract IDs from the structs
    end
    
    % Populate the existence matrix
    for i = 1:num_models
        subset_rxns = subset_models{i}.rxns;
        binaryRxnsMatrix(:, i) = ismember(ref_model.rxns, subset_rxns);
    end
    
    % Convert the existence matrix to a table
    binaryRxnsTable = array2table(binaryRxnsMatrix, 'VariableNames', model_ids);
    
    % Add the reactions list as the first column of the table
    binaryRxnsTable = addvars(binaryRxnsTable, ref_model.rxns,          'Before', 1, 'NewVariableNames', 'Reactions');
    binaryRxnsTable = addvars(binaryRxnsTable, ref_model.rxnNames, 'Before', 2, 'NewVariableNames', 'ReactionNames');
end