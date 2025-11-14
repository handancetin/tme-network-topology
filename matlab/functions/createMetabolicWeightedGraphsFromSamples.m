function createMetabolicWeightedGraphsFromSamples(model )
    % Function to create separate weighted graphs for each representative state
    
    % Step 1: Get the metabolic model
    model_name = model.id;

    modelname = strsplit(model.id, ' ');
    savename = [model.maintype + "_" + model.tissue + "_" + modelname{1}];

    fprintf('-----\n');
    fprintf('Processing model: %s\n', model_name);
    
    % Step 2: Create directed graph using S matrix
    S = full(model.S);  % Convert sparse to full matrix
    [mets, rxns] = size(S);
    
    % Initialize edge lists
    source_mets = [];
    target_mets = [];
    reaction_ids = [];
    
    % For each reaction, find substrate-product pairs
    for i = 1:rxns
        substrates = find(S(:,i) < 0);  % Negative coefficients are substrates
        products = find(S(:,i) > 0);    % Positive coefficients are products
        
        % Create edges between each substrate-product pair
        for s = 1:length(substrates)
            for p = 1:length(products)
                source_mets = [source_mets; substrates(s)];
                target_mets = [target_mets; products(p)];
                reaction_ids = [reaction_ids; i];
            end
        end
    end
    
    % Step 3: 

    randomSolutionsRHMC_filtered = load(sprintf('../output/samples/FilteredRHMC_%s.csv', savename));
    meanFlux = mean(randomSolutionsRHMC_filtered,2);

    % Create base tables
    baseNodeTable = table();
    baseNodeTable.NodeID = (1:mets)';
    baseNodeTable.MetaboliteID = model.mets;
    baseNodeTable.MetaboliteName = model.metNames;
    
    baseEdgeTable = table();
    baseEdgeTable.Source = source_mets;
    baseEdgeTable.Target = target_mets;
    baseEdgeTable.ReactionID = model.rxns(reaction_ids);
    baseEdgeTable.ReactionName = model.rxnNames(reaction_ids);
     
    % Create edge table
    edgeTable = baseEdgeTable;
    edgeTable.Flux = meanFlux(reaction_ids);
    
    % Filter edges with zero flux
    nonzero_edges = abs(edgeTable.Flux) > 1e-6;  % Using small threshold to account for numerical precision
    edgeTable = edgeTable(nonzero_edges, :);
    
    % Create adjacency matrix from remaining edges
    adj_matrix = sparse(edgeTable.Source, edgeTable.Target, ones(height(edgeTable),1), mets, mets);
    %adj_matrix = adj_matrix + adj_matrix';  % Make undirected for component analysis
    
    % Find biomass-connected component
    biomass_met_idx = find(strcmp(model.mets, 'MAM03971e'));
    if isempty(biomass_met_idx)
        warning('No biomass metabolite found in model %s', model_name); 
    end
    
    % Find connected components
    [~, ~, members] = networkComponents(adj_matrix);
    biomass_component = members == members(biomass_met_idx);
    
    % Filter nodes to biomass component
    nodeTable = baseNodeTable(biomass_component, :);
    
    % Create new consecutive node IDs
    old_to_new_id = zeros(mets, 1);
    old_to_new_id(biomass_component) = 1:sum(biomass_component);
    
    % Filter and update edges
    valid_edges = ismember(edgeTable.Source, nodeTable.NodeID) & ...
                 ismember(edgeTable.Target, nodeTable.NodeID);
    edgeTable = edgeTable(valid_edges, :);
    
    % Update edge source and target IDs
    edgeTable.Source = old_to_new_id(edgeTable.Source);
    edgeTable.Target = old_to_new_id(edgeTable.Target);
    nodeTable.NodeID = (1:height(nodeTable))';
    
    % Save state-specific tables
    output_dir = '../output/graphs/';
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    
    writetable(nodeTable, sprintf('%sNodeTable_%s.csv', output_dir, savename));
    writetable(edgeTable, sprintf('%sEdgeTable_%s.csv', output_dir, savename));
    
    fprintf('- Graph saved with %d nodes, %d edges.\n', height(nodeTable), height(edgeTable));

end

function [nComponents, sizes, members] = networkComponents(adjacency)
    % Helper function to find connected components in an undirected graph
    % Using depth-first search
    n = size(adjacency, 1);
    visited = false(n, 1);
    members = zeros(n, 1);
    nComponents = 0;
    sizes = [];
    
    for i = 1:n
        if ~visited(i)
            nComponents = nComponents + 1;
            componentSize = 0;
            stack = i;
            
            while ~isempty(stack)
                node = stack(end);
                stack(end) = [];
                
                if ~visited(node)
                    visited(node) = true;
                    members(node) = nComponents;
                    componentSize = componentSize + 1;
                    neighbors = find(adjacency(node,:));
                    stack = [stack neighbors];
                end
            end
            
            sizes(nComponents) = componentSize;
        end
    end
end