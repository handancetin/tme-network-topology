function [G] = plotReactionGraph(model)
    % Generate model name
    if isfield(model, 'maintype')
        cellname = strsplit(model.id, ' ');
        model_name = [model.maintype '_' model.tissue '_' cellname{1}];
    else
        % If model.maintype doesn't exist, use a simpler naming convention
        model_name = strrep(model.id, ' ', '_');
    end
    fprintf('Processing model: %s\n', model_name);

    % Create adjacency matrix for reactions
    S_abs = abs(model.S); % Take absolute value of stoichiometric matrix
    A = (S_abs' * S_abs) > 0; % Reactions that share metabolites
    A = A - diag(diag(A)); % Remove self-loops

    % Count connected and disconnected reactions
    total_reactions = size(A, 1);
    connected_reactions = sum(sum(A) > 0);
    disconnected_reactions = total_reactions - connected_reactions;

    % Print the count of disconnected reactions
    fprintf('Number of disconnected reactions: %d out of %d total reactions\n', ...
        disconnected_reactions, total_reactions);

    % Convert to graph
    G = graph(A, model.rxns);

    % Plot the graph
    figure;
    plot(G, 'Layout', 'force');
    title(sprintf('Reaction Graph for %s\n%d connected, %d disconnected reactions', ...
        model_name, connected_reactions, disconnected_reactions));

    % Display completion message
    fprintf('Processing complete.\n');
end