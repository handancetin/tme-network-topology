function [representative_states, indices, analysis] = fluxStatesFromRandomSamples(flux_matrix, n_states, method)
    % Extract representative states from a flux solution matrix
    %
    % Inputs:
    %   flux_matrix: Matrix of size (n_reactions x n_solutions)
    %   n_states: Number of representative states to extract (default: 5)
    %   method: Analysis method ('kmeans', 'pca', 'random', 'extreme')
    %
    % Outputs:
    %   representative_states: Matrix of size (n_reactions x n_states)
    %   indices: Indices of selected states in original matrix (for some methods)
    %   analysis: Structure containing additional analysis results
    
    if nargin < 2
        n_states = 5;
    end
    if nargin < 3
        method = 'kmeans';
    end
    
    % Transpose for clustering (solutions should be rows)
    X = flux_matrix';
    
    switch lower(method)
        case 'kmeans'
            % Use k-means clustering to find centers
            [idx, centers] = kmeans(X, n_states, 'Replicates', 5);
            representative_states = centers';
            indices = [];
            
            % Additional analysis
            analysis.cluster_sizes = accumarray(idx, 1);
            analysis.cluster_assignment = idx;
            
        case 'pca'
            % Use PCA to find principal components
            [coeff, score, latent] = pca(X);
            
            % Project data onto first n_states components
            representative_states = score(:, 1:n_states)';
            indices = [];
            
            % Additional analysis
            analysis.explained_variance = cumsum(latent) / sum(latent);
            analysis.principal_components = coeff(:, 1:n_states);
            
        case 'random'
            % Randomly select states
            indices = randperm(size(X, 1), n_states);
            representative_states = flux_matrix(:, indices);
            
            % Additional analysis
            analysis.selected_indices = indices;
            
        case 'extreme'
            % Select states with extreme flux values
            flux_norms = vecnorm(X, 2, 2);  % Calculate L2 norm for each solution
            [~, indices] = sort(flux_norms, 'descend');
            indices = indices(1:n_states);
            representative_states = flux_matrix(:, indices);
            
            % Additional analysis
            analysis.flux_norms = flux_norms(indices);
            analysis.selected_indices = indices;
            
        otherwise
            error('Unknown method specified');
    end
    
    % Common analysis for all methods
    analysis.original_mean = mean(flux_matrix, 2);
    analysis.original_std = std(flux_matrix, 0, 2);
    
    % Calculate how well each solution is represented by the states
    distances = zeros(n_states, size(flux_matrix, 2));
    for i = 1:n_states
        distances(i,:) = vecnorm(flux_matrix - representative_states(:,i), 2, 1);
    end
    analysis.min_distances = min(distances);
    analysis.mean_distance = mean(analysis.min_distances);
end