function samples = sampleForGraphs(cell_models, nStepsPerPoint, nPointsReturned)
    
    % Sample models using RHMC 
    for i = 1:length(cell_models)  
        curModel = cell_models{i};
        curModel = rmfield(curModel, {'csense'}); % Needed to fix
        cell_name = curModel.id;
        
        % Solve for maximum biomass
        sol = optimizeCbModel(curModel, 'max'); 
        
        % Constrain the cell to produce at least half of its maximum biomass
        curModel.lb(find(curModel.c)) = sol.f * 0.5;
        
        % Sample 2k points
        options.nStepsPerPoint = nStepsPerPoint;   % sampling density
        options.nPointsReturned = nPointsReturned; % number of points returned
        
        [samples, roundedPolytope, minFlux, maxFlux] = chrrSampler(curModel, nStepsPerPoint, nPointsReturned, 1, [], 0, 100);


        [~, curSamples] = sampleCbModel(curModel, [], 'CHRR', options);
        
        % Save at each iteration 
        writematrix(curSamples, fullfile('../output/samples/', sprintf('RHMC_%s.csv', cell_name))); 
    end
end
