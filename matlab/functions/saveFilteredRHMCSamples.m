function saveFilteredRHMCSamples(model)
    % Function to create separate weighted graphs for each representative state
    
    % Step 1: Get the metabolic model
    model_name = model.id;
    fprintf('-----\n');
    fprintf('Processing model: %s\n', model_name);
    
    % Step 3: Load random sampling results and normalize per 1 unit biomass
    randomSolutionsRAVEN = load(sprintf('../output/samples/RAVENCorners_%s.csv', model_name));
    randomSolutionsRAVEN = randomSolutionsRAVEN ./ randomSolutionsRAVEN(find(model.c),:); 

    randomSolutionsRHMC =load(sprintf('../output/samples/RHMC_%s.csv', model_name));
    randomSolutionsRHMC = randomSolutionsRHMC ./ randomSolutionsRHMC(find(model.c),:); 
 
    
    % Step 2: Calculate means and ranges
    rhmc_means = mean(randomSolutionsRHMC, 2);
    raven_min = min(randomSolutionsRAVEN, [], 2);
    raven_max = max(randomSolutionsRAVEN, [], 2);
    
    % Step 3: Check which means are outside range by more than 2 units
    tolerance = 5;
    below_range = rhmc_means < (raven_min - tolerance);
    above_range = rhmc_means > (raven_max + tolerance);
    means_outside_tolerance = below_range | above_range;
    reactions_outside_range = find(means_outside_tolerance);

    % 
    % % Step 3: Check which means are outside range
    % means_within_range = (rhmc_means >= raven_min) & (rhmc_means <= raven_max);
    % reactions_outside_range = find(~means_within_range);
    % 
    % Step 4: Create modified RHMC solutions with zeroed out reactions
    randomSolutionsRHMC_filtered = randomSolutionsRHMC;
    randomSolutionsRHMC_filtered(reactions_outside_range, :) = 0;



    modelname = strsplit(model.id, ' ');
    savename = [model.maintype + "_" + model.tissue + "_" + modelname{1}];
    
    % Save state-specific tables
    output_dir = '../output/samples/';
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
        
    writematrix(randomSolutionsRHMC_filtered, sprintf('%sFilteredRHMC_%s.csv', output_dir, savename)); 
    
    fprintf('- Samples saved.');

end