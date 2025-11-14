function [G] = plotMetaboliteGraph(model)

     % Generate model name
    if isfield(model, 'maintype')
        cellname = strsplit(model.id, ' ');
        model_name = [model.maintype '_' model.tissue '_' cellname{1}];
    else
        % If model.maintype doesn't exist, use a simpler naming convention
        model_name = strrep(model.id, ' ', '_');
    end
    fprintf('Processing model: %s\n', model_name);

    % Find the index of biomass in model.mets
    biomassIndex = find(strcmp(model.mets, 'MAM03971e'));
    if isempty(biomassIndex)
        error('Biomass metabolite not found in model %s.', model_name);
    end

    % Create initial adjacency matrix for metabolites
    A = (model.S ~= 0) * (model.S ~= 0)'; % Metabolites that share reactions
    A = A - diag(diag(A)); % Remove self-loops

    % Find metabolites connected to biomass
    connected = false(size(A, 1), 1);
    connected(biomassIndex) = true;
    while true
        newConnected = connected | any(A(:, connected), 2);
        if all(newConnected == connected)
            break;
        end
        connected = newConnected;
    end

    % Count disconnected metabolites
    total_metabolites = length(model.mets);
    connected_metabolites = sum(connected);
    disconnected_metabolites = total_metabolites - connected_metabolites;
    
    % Print the count of disconnected metabolites
    fprintf('Number of disconnected metabolites: %d out of %d total metabolites\n', ...
        disconnected_metabolites, total_metabolites);

    % Filter the adjacency matrix and metabolite names
    A_connected = A(connected, connected);
    metabolite_names = model.mets(connected);

    % Convert to graph
    G = graph(A_connected, metabolite_names);

    % Plot the graph
    % figure;
    % plot(G, 'Layout', 'force');
    % title(sprintf('Metabolite Graph for %s\n%d connected, %d disconnected metabolites', ...
    %     model_name, connected_metabolites, disconnected_metabolites));

    % Display completion message
    fprintf('Processing complete.\n');
end