function createMetabolicWeightedGraph(model)
    % Function to create a weighted graphs from a metabolic network model struct 
    
    % Step 1: Get the metabolic model
    model_name = model.id;
    fprintf('-----\n');
    fprintf('Processing model: %s\n', model_name);
     
    % Step 2: Create directed graph where each node is a metabolite, each
    % edge is a reaction if two metabolites are connected, using S matrix.

     % Step 2: Create directed graph using S matrix
    S = full(model.S);  % Convert sparse to full matrix
    [mets, rxns] = size(S);

    
    % Step 3: Load samples and calculate reaction statistics from random sampling
    randomSolutions = load(sprintf('../output/samples/RAVEN_%s.csv', model_name));

    % Step 4: Assign edge weights using 10 different statistical measures.

    % Step 5: Get the connected network including the biomass metabolite.

    % Step 6: Save weighted and directed graphs for each model, as
    % nodeTable and edgeTables.
         


    % Step 2: Create matrices for tracking metabolite connectivity and weights
    numMets = size(model.S, 1);
    adjacencyMatrix = zeros(numMets);
    weightMatrix = zeros(numMets);
    
    % Step 3: Process each reaction to create metabolite-to-metabolite connections
    for rxnIdx = 1:size(model.S, 2)
        if abs(fluxes(rxnIdx)) > 1e-10 % Only consider reactions with non-zero flux
            % Get substrates (negative coefficients in S matrix)
            substrates = find(model.S(:, rxnIdx) < 0);
            % Get products (positive coefficients in S matrix)
            products = find(model.S(:, rxnIdx) > 0);
            
            % Create edges between all substrates and all products
            for sub = 1:length(substrates)
                for prod = 1:length(products)
                    adjacencyMatrix(substrates(sub), products(prod)) = 1;
                    % Add flux value to weight matrix
                    weightMatrix(substrates(sub), products(prod)) = weightMatrix(substrates(sub), products(prod)) + fluxes(rxnIdx);
                end
            end
        end
    end
    
    % Create edge table with weights and reaction names
    [srcNodes, targetNodes] = find(adjacencyMatrix);  % Get all edge connections
    weights = zeros(length(srcNodes), 1);
    rxnNames = cell(length(srcNodes), 1);
    
    % Fill in edge information
    for rxnIdx = 1:size(model.S, 2)
        if abs(fluxes(rxnIdx)) > 1e-10
            substrates = find(model.S(:, rxnIdx) < 0);
            products = find(model.S(:, rxnIdx) > 0);
            
            for sub = 1:length(substrates)
                for prod = 1:length(products)
                    % Find this edge in our srcNodes/targetNodes lists
                    edgeIdx = find(srcNodes == substrates(sub) & targetNodes == products(prod));
                    if ~isempty(edgeIdx)
                        weights(edgeIdx) = weightMatrix(substrates(sub), products(prod));
                        rxnNames{edgeIdx} = model.rxns{rxnIdx};
                    end
                end
            end
        end
    end
    
    % Create edge table with numeric IDs
    edgeTable = table([srcNodes, targetNodes], weights, rxnNames, ...
        'VariableNames', {'EndNodes', 'Weight', 'EdgeName'});
    
    % Create node table with metabolite names and numeric IDs
    nodeTable = table(model.mets, (1:length(model.mets))', ...
        'VariableNames', {'Name', 'NodeID'});
    
    % Create new graph with edge and node properties
    % The EndNodes in edgeTable already contain the correct numeric indices
    G = digraph(edgeTable, nodeTable);
    
    % First convert the EndNodes column to separate source and target columns
    srcNodes = edgeTable.EndNodes(:,1);
    targetNodes = edgeTable.EndNodes(:,2);
    
    % Create a new table with separated columns
    edgeTableNew = table(srcNodes, targetNodes, edgeTable.Weight, edgeTable.EdgeName, ...
        'VariableNames', {'Source', 'Target', 'Weight', 'EdgeName'});
    
    % Save tables to CSV
    modelname = strsplit(model.id, ' ');
    savename = [model.maintype + "_" + model.tissue + "_" + modelname{1}];
    writetable(nodeTable, sprintf('../outputs/NetworkGraphs/node_table_%s.csv',savename));
    writetable(edgeTableNew, sprintf('../outputs/NetworkGraphs/edge_table_%s.csv', savename));
 
end