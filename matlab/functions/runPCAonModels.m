function [binaryRxnsTable, binaryMetsTable, binaryGenesTable] = runPCAonModels(ref_model, subset_models)

    % Extract model IDs and states 
    num_models = length(subset_models);
    model_ids = cell(1, num_models);
    model_states = cell(1, num_models); % Initialize to hold states
    for i = 1:num_models
        model_ids{i} = subset_models{i}.id; % Extract IDs from the structs
        model_states{i} = subset_models{i}.tissue; % Extract States from the structs
    end

    % Create binary matrices for rxns, mets, and genes
    [binaryRxnsMatrix, binaryRxnsTable] = createBinaryRxnsTableCFM(ref_model, subset_models);
    [binaryMetsMatrix, binaryMetsTable] = createBinaryMetsTableCFM(ref_model, subset_models);
    [binaryGenesMatrix, binaryGenesTable] = createBinaryGenesTableCFM(ref_model, subset_models);

    binaryMatrices = {binaryRxnsMatrix, binaryMetsMatrix, binaryGenesMatrix};

    % Titles and color mapping for the subplots
    titles = {'PCA of Reactions', 'PCA of Metabolites', 'PCA of Genes'};
    state_colors = struct('Normal', '#1f78b4', 'Tumor', '#e31a1c'); % Define colors for states

    % Create a figure to hold three subplots
    figure;

    % Loop through each binary matrix and generate subplots
    for i = 1:3
        % Perform PCA
        [~, score, latent] = pca(binaryMatrices{i}');
        explained = 100 * latent / sum(latent); % Calculate percentage of variance

        % Plot PCA results
        subplot(1, 3, i);
        hold on;
        for j = 1:length(model_ids)
            % Get color based on state
            color = state_colors.(model_states{j}); % Use state to determine color
            scatter(score(j, 1), score(j, 2), 150, hex2rgb(color), 'filled');
            text(score(j, 1) + 0.5, score(j, 2) + 0.01, model_ids{j}, 'FontSize', 12);
        end
        xlabel(['PC1 (' num2str(explained(1), '%.2f') '%)']);
        ylabel(['PC2 (' num2str(explained(2), '%.2f') '%)']);
        title(titles{i});
        hold off;
    end

end



function [binaryRxnsMatrix, binaryRxnsTable] = createBinaryRxnsTableCFM(ref_model, subset_models)

    num_models = length(subset_models);
    
    % Initialize the existence matrix
    binaryRxnsMatrix = zeros(length(ref_model.rxns), num_models);
    
    % Extract model names  and ensure uniqueness
    savenames = cell(1, num_models);
    for i = 1:num_models 
        model = subset_models{i};
        modelname = strsplit(model.id, ' ');
        savenames{i} = model.maintype + "_" + model.tissue + "_" + modelname{1};
    end
    savenames = [savenames{:}];


    % Populate the existence matrix
    for i = 1:num_models
        subset_rxns = subset_models{i}.rxns;
        binaryRxnsMatrix(:, i) = ismember(ref_model.rxns, subset_rxns);
    end
    
    % Collect reaction formulas with metIDs
    

    % Convert the existence matrix to a table
    binaryRxnsTable = array2table(binaryRxnsMatrix, 'VariableNames', savenames);
    
    % Add the reactions list as the first column of the table
    binaryRxnsTable = addvars(binaryRxnsTable, ref_model.rxns,          'Before', 1, 'NewVariableNames', 'Reactions');
    binaryRxnsTable = addvars(binaryRxnsTable, ref_model.rxnNames, 'Before', 2, 'NewVariableNames', 'ReactionNames');
    binaryRxnsTable = addvars(binaryRxnsTable, ref_model.subSystems, 'Before', 3, 'NewVariableNames', 'subSystems');
    binaryRxnsTable = addvars(binaryRxnsTable, ref_model.rxnFormulas, 'Before', 3, 'NewVariableNames', 'rxnFormulas');
end


function [binaryMetsMatrix, binaryMetsTable] = createBinaryMetsTableCFM(ref_model, subset_models)

    num_models = length(subset_models);
    
    % Initialize the existence matrix
    binaryMetsMatrix = zeros(length(ref_model.mets), num_models);
    
     % Extract model names  and ensure uniqueness
    savenames = cell(1, num_models);
    for i = 1:num_models 
        model = subset_models{i};
        modelname = strsplit(model.id, ' ');
        savenames{i} = [model.maintype + "_" + model.tissue + "_" + modelname{1}];
    end
    savenames = [savenames{:}];
    
    % Populate the existence matrix
    for i = 1:num_models
        subset_mets = subset_models{i}.mets;
        binaryMetsMatrix(:, i) = ismember(ref_model.mets, subset_mets);
    end
    
    % Convert the existence matrix to a table
    binaryMetsTable = array2table(binaryMetsMatrix, 'VariableNames', savenames);
    
    % Add the reactions list as the first column of the table
    binaryMetsTable = addvars(binaryMetsTable, ref_model.mets,          'Before', 1, 'NewVariableNames', 'Metabolites');
    binaryMetsTable = addvars(binaryMetsTable, ref_model.metNames, 'Before', 2, 'NewVariableNames', 'MetaboliteNames');
end


function [binaryGenesMatrix, binaryGenesTable] = createBinaryGenesTableCFM(ref_model, subset_models)

    num_models = length(subset_models);
    
    % Initialize the existence matrix
    binaryGenesMatrix = zeros(length(ref_model.genes), num_models);
    
     % Extract model names  and ensure uniqueness
    savenames = cell(1, num_models);
    for i = 1:num_models 
        model = subset_models{i};
        modelname = strsplit(model.id, ' ');
        savenames{i} = [model.maintype + "_" + model.tissue + "_" + modelname{1}];
    end
    savenames = [savenames{:}];
    
    % Populate the existence matrix
    for i = 1:num_models
        subset_genes = subset_models{i}.genes;
        binaryGenesMatrix(:, i) = ismember(ref_model.genes, subset_genes);
    end
    
    % Convert the existence matrix to a table
    binaryGenesTable = array2table(binaryGenesMatrix, 'VariableNames', savenames);
    
    % Add the reactions list as the first column of the table
    binaryGenesTable = addvars(binaryGenesTable, ref_model.genes,          'Before', 1, 'NewVariableNames', 'Genes');
    binaryGenesTable = addvars(binaryGenesTable, ref_model.geneShortNames, 'Before', 2, 'NewVariableNames', 'GeneShortNames');
   
    % Loop through each gene
    usedReactions = cell(length(ref_model.genes), 1);
    for i = 1:length(ref_model.genes)
        % Get the current gene ID
        currentGene = ref_model.genes{i};
        
        % Find reactions that contain this gene in their grRules
        % We need to search for the exact gene ID (with word boundaries)
        rxnIndices = [];
        for j = 1:length(ref_model.grRules)
            % Check if the current gene is in this rule
            if contains(ref_model.grRules{j}, currentGene)
                rxnIndices = [rxnIndices; j];
            end
        end
        
        % If there are reactions using this gene
        if ~isempty(rxnIndices)
            % Get the reaction IDs
            geneRxns = ref_model.rxns(rxnIndices);
            
            % Join them with commas
            usedReactions{i} = strjoin(geneRxns, ', ');
        else
            % If no reactions use this gene, store an empty string
            usedReactions{i} = '';
        end
    end
    binaryGenesTable = addvars(binaryGenesTable, usedReactions, 'Before', 3, 'NewVariableNames', 'UsedReactions');
end