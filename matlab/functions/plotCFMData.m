
function plotCFMData(expr_data)


tmp = quantile(expr_data.cpm_scaled, 0.90);
expr_data.cpm_scaled(expr_data.cpm_scaled <tmp) = 0;


% Define a struct array mapping each maintype to its corresponding color
maintype_colors = struct('type', {'Cancer', 'Fibroblast', 'Macrophage'}, ...
                                       'color', {'#d7191c', '#bababa', '#2c7bb6'});

% Initialize the colors matrix
colors = zeros(length(expr_data.maintypes), 3);

% Loop through each maintype and assign colors
for i = 1:length(maintype_colors)
    mask = strcmp(expr_data.maintypes, maintype_colors(i).type);
    colors(mask, :) = repmat(hex2rgb(maintype_colors(i).color), sum(mask), 1);
end

% Plot each box with its own color
figure()
hold on; % Keep all boxes on the same plot
numBoxes = size(expr_data.cpm_scaled, 2); % Number of boxes
boxHandles = gobjects(1, numBoxes); % Array to store box handles
maxY = 0;
for i = 1:numBoxes
    % Extract data for the current box
    boxData = log2(expr_data.cpm_scaled(:, i)+1);
    
    % Create the boxchart for the current box
    boxHandles(i) = boxchart(i * ones(size(boxData)), boxData, ...
                             'JitterOutliers', 'on', ...
                             'MarkerStyle', '.', ...
                             'MarkerSize', 15, ...
                             'MarkerColor', colors(i, :), ...
                             'BoxFaceColor', colors(i, :), ...
                             'BoxFaceAlpha', .8, ...
                             'BoxEdgeColor', 'black');
    maxY = max(maxY, max(boxData));
end
ylabel('log_2(cpm+1)');
xticks(1:numBoxes); % Set XTick positions
xlim([0, numBoxes+1])
ylim([-0.5, maxY+0.5])

xticklabels(expr_data.clusters);

hold off;



% Perform PCA
[~, score, latent] = pca(expr_data.cpm_scaled');
explained = 100 * latent / sum(latent); % Calculate percentage of variance

figure()
hold on;
for j = 1:length(expr_data.clusters)
    scatter(score(j,1), score(j,2), 150, 'k', 'filled');
    text(score(j,1) + 0.5, score(j,2) + 0.01, expr_data.clusters{j}, 'FontSize', 12);
end
xlabel(['PC1 (' num2str(explained(1), '%.2f') '%)']);
ylabel(['PC2 (' num2str(explained(2), '%.2f') '%)']); 
hold off;

end