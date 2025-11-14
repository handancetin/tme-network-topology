
function [binaryGenesMatrix, binaryGenesTable] = createBinaryGenesTableCFM(ref_model, subset_models)

    num_models = length(subset_models);
    
    % Initialize the existence matrix
    binaryGenesMatrix = zeros(length(ref_model.genes), num_models);
    
    % Extract model IDs and ensure uniqueness
    model_ids = cell(1, num_models);
    for i = 1:num_models
        model_ids{i} = subset_models{i}.id; % Extract IDs from the structs
    end
    
    % Populate the existence matrix
    for i = 1:num_models
        subset_genes = subset_models{i}.genes;
        binaryGenesMatrix(:, i) = ismember(ref_model.genes, subset_genes);
    end
    
    % Convert the existence matrix to a table
    binaryGenesTable = array2table(binaryGenesMatrix, 'VariableNames', model_ids);
    
    % Add the reactions list as the first column of the table
    binaryGenesTable = addvars(binaryGenesTable, ref_model.genes,          'Before', 1, 'NewVariableNames', 'Genes');
    binaryGenesTable = addvars(binaryGenesTable, ref_model.geneShortNames, 'Before', 2, 'NewVariableNames', 'GeneShortNames');
end