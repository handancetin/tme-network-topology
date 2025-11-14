function [solutions, goodRxns] = sampleUsingRavenSampler(cell_models, nSamples)
    
    % Sample models using RHMC 
    for i = 1:length(cell_models)  
        curModel = cell_models{i};
        cell_name = curModel.id;
        biomass_index = find(curModel.c);
        
        % Solve for maximum biomass
        curModel = rmfield(curModel, {'csense'}); % Needed to fix
        sol = solveLP(curModel);
        
        % Constrain the cell to produce at least half of its maximum biomass
        curModel.lb(biomass_index) = sol.f * 0.25;
        
        % Sample points
        [solutions, goodRxns] = randomSampling(curModel,nSamples,false,false,false,[],true);

        % Normalize samples
        %solutions = solutions ./ solutions(biomass_index,:);

        % Save at each iteration 
        writematrix(solutions, fullfile('../output/samples/', sprintf('RAVENCorners_%s.csv', cell_name))); 
    end
end
