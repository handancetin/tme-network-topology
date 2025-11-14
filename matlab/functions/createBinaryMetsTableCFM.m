
function [binaryMetsMatrix, binaryMetsTable] = createBinaryMetsTableCFM(ref_model, subset_models)

    num_models = length(subset_models);
    
    % Initialize the existence matrix
    binaryMetsMatrix = zeros(length(ref_model.mets), num_models);
    
    % Extract model IDs and ensure uniqueness
    model_ids = cell(1, num_models);
    for i = 1:num_models
        model_ids{i} = subset_models{i}.id; % Extract IDs from the structs
    end
    
    % Populate the existence matrix
    for i = 1:num_models
        subset_mets = subset_models{i}.mets;
        binaryMetsMatrix(:, i) = ismember(ref_model.mets, subset_mets);
    end
    
    % Convert the existence matrix to a table
    binaryMetsTable = array2table(binaryMetsMatrix, 'VariableNames', model_ids);
    
    % Add the reactions list as the first column of the table
    binaryMetsTable = addvars(binaryMetsTable, ref_model.mets,          'Before', 1, 'NewVariableNames', 'Metabolites');
    binaryMetsTable = addvars(binaryMetsTable, ref_model.metNames, 'Before', 2, 'NewVariableNames', 'MetaboliteNames');
end