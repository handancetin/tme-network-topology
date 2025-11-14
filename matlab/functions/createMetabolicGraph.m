function adjacencyMatrix = createMetabolicGraph(model, plotGraph, saveGraph)
    % Function to create a graph from a metabolic network model struct.
    %
    % Input:
    %   model - a struct representing the metabolic model with fields:
    %       .id - model ID
    %       .maintype - model maintype
    %       .tissue - model tissue type
    %       .mets - cell array of metabolite names
    %       .S - stoichiometric matrix
    %
    % Output:
    %   adjSim - adjacency matrix of metabolites connected in the network
    
    % Generate model name
    cellname = strsplit(model.id, ' ');
    model_name = [model.maintype '_' model.tissue '_' cellname{1}];
    
    fprintf('Processing model: %s\n', model_name);
    
    try
        % Step 1: Find the index of biomass in model.mets
        biomassIndex = find(strcmp(model.mets, 'MAM03971e'));
        if isempty(biomassIndex)
            warning('Biomass metabolite not found in model %s. Skipping.', model_name);
            adjSim = [];
            return;
        end

        % Step 2: Create initial adjacency matrix for metabolites
        A = (model.S ~= 0) * (model.S ~= 0)';  % Metabolites that share reactions
        A = A - diag(diag(A));  % Remove self-loops

        % Step 3: Find metabolites connected to biomass
        connected = false(size(A, 1), 1);
        connected(biomassIndex) = true;
        while true
            newConnected = connected | any(A(:, connected), 2);
            if all(newConnected == connected)
                break;
            end
            connected = newConnected;
        end

        % Step 4: Filter the adjacency matrix and metabolite names
        A_connected = A(connected, connected);
        metabolite_names = model.mets(connected);

        % Step 5: Convert to graph
        G = graph(A_connected, metabolite_names);

       % Plot the graph
       if plotGraph == 1
            figure;
            %p = plot(G, 'NodeLabel', metabolite_names);
            p = plot(G);
            
            % Customize appearance
            p.MarkerSize = 6;       % Set node size
            p.LineWidth = 1;      % Set edge thickness
            
            % Get the number of nodes and edges
            numNodes = numnodes(G);
            numEdges = numedges(G);
            
            % Add title with number of nodes and edges
            title(sprintf('Metabolic Network Graph: %d Nodes, %d Edges', numNodes, numEdges));
       end

        % Step 6: Get the adjacency matrix from the graph
        adjacencyMatrix = full(adjacency(G));
          
        if saveGraph == 1
            save_filename = fullfile('../outputs/graphs', ['graph_' model_name '.mat']);
            save(save_filename, 'adjacencyMatrix'); 
            fprintf('Saved graph for %s\n', model_name);


            % Save metabolite names to 'nodes_' model_name '.txt'
            nodes_filename = fullfile('../outputs/graphs', ['nodes_' model_name '.txt']);
            fileID = fopen(nodes_filename, 'w');
        
            % Write each metabolite name to a new line
            for i = 1:length(metabolite_names)
                fprintf(fileID, '%s\n', metabolite_names{i});
            end
        
            % Close the file
            fclose(fileID);
            fprintf('Saved nodes for %s\n', model_name);

        end

    catch ME
        warning('Error processing model %s: %s', model_name, ME.message);
        adjacencyMatrix = [];
    end
end
