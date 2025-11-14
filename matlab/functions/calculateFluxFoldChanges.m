function results = calculateFluxFoldChanges(dataTable, modelName1, modelName2, varargin)

    % Only investigate common rows, remove if one of the flux is 0
    dataTable((dataTable{:,3} == 0) | (dataTable{:,4} == 0), :) = [];
    flux1 = dataTable.(modelName1);
    flux2 = dataTable.(modelName2);
    
    % Calculate fold change with handling for zero and negative values
    epsilon = 1e-10; % Small value to avoid division by zero
    fold_change = (abs(flux2) + epsilon) ./ (abs(flux1) + epsilon);
    
    % Adjust fold change sign based on the original flux signs
    fold_change = fold_change .* sign(flux2 .* flux1);
    
    % Calculate log2 fold change
    log2_fold_change = log2(abs(fold_change)) .* sign(fold_change);
    
    % Create results table
    results = table((1:length(flux1))', dataTable.Reactions, dataTable.ReactionNames, flux1, flux2, fold_change, log2_fold_change, ...
                             'VariableNames', {'rxnID', 'Reactions', 'ReactionNames', 'Flux1', 'Flux2', 'FoldChange', 'Log2FoldChange'});
    
    % Add flags for zero and negative flux
    results.ZeroFlux1 = abs(flux1) < epsilon;
    results.ZeroFlux2 = abs(flux2) < epsilon;
    results.NegativeFlux1 = flux1 < 0;
    results.NegativeFlux2 = flux2 < 0;
    
    % Sort results by absolute log2 fold change
    results.AbsLog2FoldChange = abs(results.Log2FoldChange);
    results = sortrows(results, 'AbsLog2FoldChange', 'descend');
    
    % Display top 10 rows
    disp('Top 10 rows by absolute log2 fold change:');
    disp(results(1:10,:)); 
    
    % Filter for significant changes based on fold change
    fc_threshold = 2; % fold change threshold 
    significant_changes = results(results.AbsLog2FoldChange > log2(fc_threshold), :);
    disp('Number of significant changes:');
    disp(height(significant_changes));
    
    
    % Create a scatter plot with log2 fold changes
    figure;
                scatter(results.Flux1, results.Flux2, 40, results.Log2FoldChange, 'filled');
                hold on;
    
                % Calculate axis limits
                minFlux = min(min(results.Flux1), min(results.Flux2));
                maxFlux = max(max(results.Flux1), max(results.Flux2));
                
                % Plot the diagonal line
                plot([minFlux, maxFlux], [minFlux, maxFlux], 'r--');
                
                % Label points with |log2 fold change| > 15
                high_change_idx = abs(results.Log2FoldChange) > 15;
                text(results.Flux1(high_change_idx), results.Flux2(high_change_idx), ...
                    results.Reactions(high_change_idx), 'FontSize', 12, 'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'bottom', 'Rotation', 0);
                
                % Set figure details
                xlabel('Flux in Model 1');
                ylabel('Flux in Model 2');
                title('Reaction Flux Comparison with Log2 Fold Change (Names are shown for FC>15)');
                colorbar;
                hold off;

end