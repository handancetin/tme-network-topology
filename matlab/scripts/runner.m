initCustomSettings()

% generated results will be saved in this directiory
directory_output = '../output';

% Load the expression data
exprData = loadCFMExprData('/Users/handan/Documents/CSBL/project-files/celltype-gem-interactions/r/output/averaged_lognorm_expressions_cpm.csv'); 
plotCFMData(exprData)

%% run ftINIT
load('prepData_essentialTasks.mat')

CFM_models = cell(length(exprData.clusters), 1);
ftinit_details = cell(length(exprData.clusters), 1);
for i = 1:length(CFM_models)
     [CFM_models{i}, ftinit_details{i}]  = createClusterModels(exprData.genes, exprData.clusters(i), exprData.maintypes(i), exprData.tissues(i), exprData.cpm_scaled(:, i), 0.50, prepData_essentialTasks);
end
save('CFM_models.mat', 'CFM_models')
save('ftinit_details.mat', 'CFM_models')


%% Delete
%% run ftINIT
load('prepData_essentialTasks.mat')

CFM_models_90 = cell(length(exprData.clusters), 1);
ftinit_details_90 = cell(length(exprData.clusters), 1);
for i = 1:length(CFM_models)
     [CFM_models_90{i}, ftinit_details_90{i}]  = createClusterModels(exprData.genes, exprData.clusters(i), exprData.maintypes(i), exprData.tissues(i), exprData.cpm_scaled(:, i), 0.90, prepData_essentialTasks);
end
save('CFM_models_90.mat', 'CFM_models_90')
save('ftinit_details_90.mat', 'CFM_models_90')

%% save rxn Tables
load CFM_models
load hg_rs_model.mat

% Add reaction formulas into struct
hg_rs_model.rxnFormulas = printRxnFormula(hg_rs_model, 'rxnAbbrList', hg_rs_model.rxns, 'printFlag', 0);
  
[binaryRxnsTable,binaryMetsTable, binaryGenesTable] = runPCAonModels(hg_rs_model, CFM_models);

% Save each table as CSV
writetable(binaryRxnsTable, fullfile(directory_output, 'binaryRxnsTable.csv'));
writetable(binaryMetsTable, fullfile(directory_output, 'binaryMetsTable.csv'));
writetable(binaryGenesTable, fullfile(directory_output, 'binaryGenesTable.csv'));

%% Sample Models
sampleUsingRavenSampler(CFM_models, 1000); % done
sampleUsingRHMC(CFM_models, 200, 2000) % done

% Save RAVEN filtered RHMC samples future analysis
for i = 1:length(CFM_models)
  saveFilteredRHMCSamples(CFM_models{i});
end 
   
% Save metabolite graphs 
for i = 1:length(CFM_models)
  createMetabolicWeightedGraphsFromSamples(CFM_models{i});
end

 
%% Optimize for growth
load CFM_models
load hg_rs_model.mat

binaryRxnsTable = readtable(fullfile(directory_output, 'binaryRxnsTable.csv'));
[maxGrowths, solutionsTable]  = runMinFBAonModels(CFM_models, binaryRxnsTable, 1);
writetable(solutionsTable, '../output/simulations/FBAsolutionsTable_solveLP_minFlux1.csv'); 


%% Run single reaction deletion for all

[knockoutGrowth, knockoutRatios] = runReactionDeletiononModels(CFM_models, binaryRxnsTable, 'FBA');
writetable(knockoutGrowth, '../output/simulations/singleReaction_knockoutGrowth_FBA.csv'); 
writetable(knockoutRatios, '../output/simulations/singleReaction_knockoutRatios_FBA.csv'); 
 
[knockoutGrowth_MOMA, knockoutRatios_MOMA] = runReactionDeletiononModels(CFM_models, binaryRxnsTable, 'MOMA');
writetable(knockoutGrowth_MOMA, '../output/simulations/singleReaction_knockoutGrowth_MOMA.csv'); 
writetable(knockoutRatios_MOMA, '../output/simulations/singleReaction_knockoutRatios_MOMA.csv'); 


%% Run single gene deletion for all

[knockoutGeneGrowth, knockoutGeneRatios] = runGeneDeletionOnModels(CFM_models, binaryGenesTable, 'FBA');
writetable(knockoutGeneGrowth, '../output/simulations/singleGene_knockoutGrowth_FBA.csv'); 
writetable(knockoutGeneRatios, '../output/simulations/singleGene_knockoutRatios_FBA.csv'); 
 
[knockoutGeneGrowth, knockoutGeneRatios] = runGeneDeletionOnModels(CFM_models, binaryGenesTable, 'MOMA');
writetable(knockoutGeneGrowth, '../output/simulations/singleGene_knockoutGrowth_MOMA.csv'); 
writetable(knockoutGeneRatios, '../output/simulations/singleGene_knockoutRatios_MOMA.csv');

%% Try knockouting leukotriene and arachidonic acid metabolism	in MARCO

for i= 1:2:length(CFM_models)
    normal_model = CFM_models{i, 1}  ;
    tumor_model = CFM_models{i+1, 1}  ;
    
    normal_model_ko = normal_model;
    tumor_model_ko = tumor_model;
    
    % tumor_model_ko.lb(find(contains(tumor_model_ko.subSystems, 'Porphyrin metabolism'))) = 0;
    % tumor_model_ko.ub(find(contains(tumor_model_ko.subSystems, 'Porphyrin metabolism'))) = 0;

    tumor_model_ko.lb(find(contains(tumor_model_ko.rxns, 'MAR02254'))) = 0;
    tumor_model_ko.ub(find(contains(tumor_model_ko.rxns, 'MAR02254'))) = 0;
    
    % tumor_model_ko.lb(find(contains(tumor_model_ko.subSystems, 'Arachidonic acid metabolism'))) = 0;
    % tumor_model_ko.ub(find(contains(tumor_model_ko.subSystems, 'Arachidonic acid metabolism'))) = 0;
     
    % normal_model_ko.lb(find(contains(normal_model_ko.subSystems, 'Porphyrin metabolism'))) = 0;
    % normal_model_ko.ub(find(contains(normal_model_ko.subSystems, 'Porphyrin metabolism'))) = 0;


    normal_model_ko.lb(find(contains(normal_model_ko.rxns, 'MAR02254'))) = 0;
    normal_model_ko.ub(find(contains(normal_model_ko.rxns, 'MAR02254'))) = 0;
    
    % normal_model_ko.lb(find(contains(normal_model_ko.subSystems, 'Arachidonic acid metabolism'))) = 0;
    % normal_model_ko.ub(find(contains(normal_model_ko.subSystems, 'Arachidonic acid metabolism'))) = 0;
      
    fprintf('%s\n', normal_model.id)
    fprintf('-- WT model, in Normal tissue: %f \n', solveLP(normal_model, 1).f)
    fprintf('-- KO model, in Normal tissue: %f \n', solveLP(normal_model_ko, 1).f)
    fprintf('-- WT model, in Tumor tissue: %f \n', solveLP(tumor_model, 1).f)
    fprintf('-- KO model, in Tumor tissue: %f \n', solveLP(tumor_model_ko, 1).f)
    fprintf('-- Differences: %f and  %f \n', solveLP(normal_model, 1).f - solveLP(normal_model_ko, 1).f, solveLP(tumor_model, 1).f - solveLP(tumor_model_ko, 1).f)
end


sol_wt = solveLP(marco_tumor, 1).x
sol_ko = solveLP(tumor_model_ko, 1).x

sol_ratio = [sol_wt, sol_ko, abs(sol_ko./(sol_wt+0.000001))];


solutionsTable(find(contains(solutionsTable.subSystems, 'Arachidonic acid metabolism')), :)