%% Set the working environment

% NOTE: 
% make sure toolboxes are in the MATLAB path: RAVEN, COBRA, GECKO
% then, set solvers of toolboxes to gurobi (academically free)

% generated results will be saved in this directiory
directory_output = '../output';

% personal preference for figures, can be skipped
set(0, 'DefaultFigureWindowStyle','docked')
set(0,'defaultFigureColor', [1 1 1])
set(groot,'defaultAxesFontName','Arial')
set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(groot,'defaultAxesGridAlpha', 0.05)
set(groot,'defaultAxesGridColor', [0 0 0])
set(groot,'defaultAxesColor', [1,1,1])%[0.99 0.99 0.99])
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesTickDirMode', 'manual');
set(groot,'defaultAxesBox', 'on')
set(groot,'defaultAxesLineWidth', 1)


%% Read pseudobulk expression profiles

% NOTE:
% n < 50 models can be removed in this section for compuatation speed. 
% simulation results of those models will not be used in statistical 
% analyses, however, this script generates models for those, too.

% load expression profiles (obtained with r)
expressions = readtable('../../r/output/averaged_lognorm_expressions_cpm.csv');

% fix column names
name_mapping = {
    'Var1', 'Genes';
    'Malignant.cells_N', 'Malignant cells (N)';
    'Malignant.cells_T', 'Malignant cells (T)';
    'FGFR2..fibroblasts_N', 'FGFR2+ fibroblasts (N)';
    'FGFR2..fibroblasts_T', 'FGFR2+ fibroblasts (T)';
    'ICAM1..telocytes_N', 'ICAM1+ telocytes (N)';
    'ICAM1..telocytes_T', 'ICAM1+ telocytes (T)';
    'FAP..fibroblasts_N', 'FAP+ fibroblasts (N)';
    'FAP..fibroblasts_T', 'FAP+ fibroblasts (T)';
    'ICAM1..telocytes_N.1', 'ICAM1- telocytes (N)';
    'ICAM1..telocytes_T.1', 'ICAM1- telocytes (T)';
    'MFAP5..myofibroblasts_N', 'MFAP5+ myofibroblasts (N)';
    'MFAP5..myofibroblasts_T', 'MFAP5+ myofibroblasts (T)';  % n < 50
    'CD73..fibroblasts_N', 'CD73+ fibroblasts (N)';
    'CD73..fibroblasts_T', 'CD73+ fibroblasts (T)';
    'DES..myofibroblasts_N', 'DES+ myofibroblasts (N)';
    'DES..myofibroblasts_T', 'DES+ myofibroblasts (T)';
    'THBS1..Macrophage_N', 'THBS1+ Macrophage (N)';
    'THBS1..Macrophage_T', 'THBS1+ Macrophage (T)';
    'VCAN..Monocyte_N', 'VCAN+ Monocyte (N)';
    'VCAN..Monocyte_T', 'VCAN+ Monocyte (T)';
    'MARCO..Macrophage_N', 'MARCO+ Macrophage (N)';  % n < 50
    'MARCO..Macrophage_T', 'MARCO+ Macrophage (T)';
};
expressions.Properties.VariableNames = name_mapping(:, 2); clear name_mapping;

% create struct for context-specific model reconstruction
expr_struct.genes = expressions.Genes; % gene names in ENSG format
expr_struct.tissues = expressions.Properties.VariableNames(2:end); % names
expr_struct.levels = table2array(expressions(:, 2:end)); % cpm values
expr_struct.maintype = [repmat("Cancer", 1, 2),  ...
                        repmat("Fibroblast", 1, 14), ...
                        repmat("Macrophage", 1, 6)]'; % main type

% save number of samples (this variable will be used in loops)
numSamp = length(expr_struct.tissues);

%% Load generic human model Human-GEM & prepare for ftINIT
%changeCobraSolver('gurobi', 'all');
setRavenSolver('gurobi');

% Load Human-GEM
directory_humanGEM = '/Users/hcetin/Documents/CSBL/github-projects/third-party/Human-GEM-1.18.0';
load(sprintf('%s/model/Human-GEM.mat', directory_humanGEM));
model = ihuman;

% parse tasks
tasks_Essential = parseTaskList(sprintf('%s/data/metabolicTasks/metabolicTasks_Essential.txt', directory_humanGEM));
tasks_Cellfie = parseTaskList(sprintf('%s/data/metabolicTasks/metabolicTasks_CellfieConsensus.txt', directory_humanGEM));

% cancer (malignant) cells
tasksCancer  = tasks_Essential;
tasksCancer(end+1) = tasks_Cellfie(5);   % 'ATP regeneration from glucose (normoxic conditions) - glycolysis + krebs cycle'
tasksCancer(end+1) = tasks_Cellfie(6);   % 'ATP generation from glucose (hypoxic conditions) - glycolysis'
tasksCancer(end+1) = tasks_Cellfie(7);   % 'Reactive oxygen species detoxification (H2O2 to H2O)'
tasksCancer(end+1) = tasks_Cellfie(32);  % 'Glucose to lactate conversion'
tasksCancer(end+1) = tasks_Cellfie(88);  % 'Glutamine synthesis'
tasksCancer(end+1) = tasks_Cellfie(89);  % 'Glutamine degradation'
tasksCancer(end+1) = tasks_Cellfie(90);  % 'Glutaminolysis (glutamine to lactate)'
tasksCancer(end+1) = tasks_Cellfie(111); % 'Serine synthesis'
tasksCancer(end+1) = tasks_Cellfie(136); % 'Cholesterol synthesis'
tasksCancer(end+1) = tasks_Cellfie(149); % 'Palmitate synthesis'
tasksCancer(end+1) = tasks_Cellfie(150); % 'Palmitate degradation'
tasksCancer(end+1) = tasks_Cellfie(161); % 'Arachidonate synthesis'
tasksCancer(end+1) = tasks_Cellfie(169); % 'Synthesis of thromboxane from arachidonate'

prepData_cancer = prepHumanModelForftINIT(model, false, tasksCancer, sprintf('%s/model/reactions.tsv', directory_humanGEM));
save('checkpoints/prepData_cancer.mat', 'prepData_cancer')

% fibroblasts
tasksFibroblast = tasks_Essential;
tasksFibroblast(end+1) = tasks_Cellfie(5);   % 'ATP regeneration from glucose (normoxic conditions) - glycolysis + krebs cycle'
tasksFibroblast(end+1) = tasks_Cellfie(7);   % 'Reactive oxygen species detoxification (H2O2 to H2O)'
tasksFibroblast(end+1) = tasks_Cellfie(32);  % 'Glucose to lactate conversion'
tasksFibroblast(end+1) = tasks_Cellfie(88);  % 'Glutamine synthesis'
tasksFibroblast(end+1) = tasks_Cellfie(89);  % 'Glutamine degradation'
tasksFibroblast(end+1) = tasks_Cellfie(117); % 'Proline synthesis'
tasksFibroblast(end+1) = tasks_Cellfie(118); % 'Proline degradation'
tasksFibroblast(end+1) = tasks_Cellfie(161); % 'Arachidonate synthesis'

prepData_fibroblast = prepHumanModelForftINIT(model, false, tasksFibroblast, sprintf('%s/model/reactions.tsv', directory_humanGEM));
save('checkpoints/prepData_fibroblast.mat', 'prepData_fibroblast')

% macrophages
tasksMacrophage = tasks_Essential;
tasksMacrophage(end+1) = tasks_Cellfie(5);   % 'ATP regeneration from glucose (normoxic conditions) - glycolysis + krebs cycle'
tasksMacrophage(end+1) = tasks_Cellfie(7);   % 'Reactive oxygen species detoxification (H2O2 to H2O)'
tasksMacrophage(end+1) = tasks_Cellfie(62);  % 'Arginine synthesis'
tasksMacrophage(end+1) = tasks_Cellfie(63);  % 'Arginine degradation'
tasksMacrophage(end+1) = tasks_Cellfie(88);  % 'Glutamine synthesis'
tasksMacrophage(end+1) = tasks_Cellfie(89);  % 'Glutamine degradation'
tasksMacrophage(end+1) = tasks_Cellfie(121); % 'Tryptophan degradation'
tasksMacrophage(end+1) = tasks_Cellfie(123); % 'Synthesis of kynate from tryptophan'
tasksMacrophage(end+1) = tasks_Cellfie(124); % 'Synthesis of L-kynurenine from tryptophan'
tasksMacrophage(end+1) = tasks_Cellfie(161); % 'Arachidonate synthesis'

 = prepHumanModelForftINIT(model, false, tasksMacrophage, sprintf('%s/model/reactions.tsv', directory_humanGEM));
save('checkpoints/prepData_macrophage.mat', 'prepData_macrophage')


%% run ftINIT to obtain models
cfm_models = cell(numSamp, 1);
ftinit_info = struct();
for i = 1:numSamp 
    disp(['Model: ' num2str(i) ' of ' num2str(numSamp)])
    if expr_struct.maintype{i} == "Cancer"
        prepData = prepData_cancer;
    elseif expr_struct.maintype{i} == "Fibroblast"
        prepData = prepData_fibroblast;
    elseif expr_struct.maintype{i} == "Macrophage"
        prepData = prepData_macrophage;
    end
    % run ftINIT
    [cfm_models{i}, ~, ...
    ftinit_info.addedRxnsForTasks{i}, ...
    ftinit_info.deletedRxnsInINIT{i}, ...
    ftinit_info.fullMipRes{i}] = ftINIT(prepData, ...
    expr_struct.tissues{i}, [], [], expr_struct, [], getHumanGEMINITSteps('1+0'), false, true, [], false);
    cfm_models{i}.id = expr_struct.tissues{i};
end
%save('checkpoints/cfm_models.mat', 'cfm_models')
%save('checkpoints/ftinit_info.mat', 'ftinit_info')
load('checkpoints/cfm_models.mat')

%% testing metabolic tasks for generated models
% test essential cellular tasks - this is just to double check if ftINIT worked fine
taskReportsEssential = cell(numSamp, 1);
taskOkEssential = zeros(length(tasks_Essential), numSamp);
for i = 1:numSamp 
    disp(['Model: ' num2str(i) ' of ' num2str(numSamp)])
        taskReportsEssential{i} = checkTasks(closeModel(cfm_models{i}),[],false,false,false,tasks_Essential); 
        taskOkEssential(:,i) = taskReportsEssential{i}.ok;
end
saveTable = array2table(taskOkEssential, 'VariableNames', expr_struct.tissues);
saveTable = [table({tasks_Essential.description}.', 'VariableNames', {'Task'}), saveTable];
writetable(saveTable, sprintf('%s/metabolic-tasks/essential.csv', directory_output), 'WriteRowNames', false);

% test celfie tasks
taskReportsCellfie = cell(numSamp, 1);
taskOkCellfie = zeros(length(tasks_Cellfie), numSamp);
for i = 1:numSamp 
    disp(['Model: ' num2str(i) ' of ' num2str(numSamp)])
        taskReportsCellfie{i} = checkTasks(closeModel(cfm_models{i}),[],false,false,false,tasks_Cellfie); 
        taskOkCellfie(:,i) = taskReportsCellfie{i}.ok;
end
saveTable = array2table(taskOkCellfie, 'VariableNames', expr_struct.tissues);
saveTable = [table({tasks_Cellfie.description}.', 'VariableNames', {'Task'}), saveTable];
writetable(saveTable, sprintf('%s/metabolic-tasks/cellfie.csv', directory_output), 'WriteRowNames', false);


% test celfie tasks for generic model: for comparison
taskReportsCellfie_ihuman = checkTasks(closeModel(ihuman),[],false,false,false,tasks_Cellfie); 
taskOkCellfie_ihuman      = taskReportsCellfie_ihuman.ok;
saveTable = array2table(taskOkCellfie_ihuman, 'VariableNames', {'Human-Gem'});
saveTable = [table({tasks_Cellfie.description}.', 'VariableNames', {'Task'}), saveTable];
writetable(saveTable, sprintf('%s/metabolic-tasks/cellfie_ihuman.csv', directory_output), 'WriteRowNames', false);

%% constraint by media and reduce models

% Ham's media composition
availableMets = {'glucose'
                 'arginine'
                 'histidine'
                 'lysine'
                 'methionine'
                 'phenylalanine'
                 'tryptophan'
                 'tyrosine'
                 'alanine'
                 'glycine'
                 'serine'
                 'threonine'
                 'aspartate'
                 'glutamate'
                 'asparagine'
                 'glutamine'
                 'isoleucine'
                 'leucine'
                 'proline'
                 'valine'
                 'cysteine'
                 'thiamin'
                 'hypoxanthine'
                 'folate'
                 'biotin'
                 'pantothenate'
                 'choline'
                 'inositol'
                 'nicotinamide'
                 'pyridoxine'
                 'riboflavin'
                 'thymidine'
                 'aquacob(III)alamin'
                 'lipoic acid'
                 'sulfate'
                 'linoleate'
                 'linolenate'
                 'O2'
                 'H2O'
                 'retinoate'
                 'Fe2+'
                 'Pi'
                 'alpha-tocopherol'
                 'gamma-tocopherol'};

cfm_models_constrained = cell(numSamp, 1);
fva_Active = cell(numSamp, 1);
fva_Max = cell(numSamp, 1);
fva_Min = cell(numSamp, 1);
for i = 1:numSamp  
    disp(['Constraining model: ' num2str(i) ' of ' num2str(numSamp)])
    curModel = cfm_models{i}; 
    % Only allow select metabolites
    curModel.lb(find(findExcRxns(curModel, 0, 0))) = 0; % turn off all uptakes
    curModel.lb(fun_findUptakeRxnsFromMets(curModel, availableMets)) = -1000;  
    % run FVA to reduce models and enforce bounds
    [cfm_models_constrained{i}, fva_Active{i}, fva_Max{i}, fva_Min{i}] = reduceModel(curModel, 1e-6, 0, 0, 0, 1, 1);
end

% save constrained models
save('checkpoints/cfm_models_constrained.mat', 'cfm_models_constrained')

% save fva results
fva_results.names = expr_struct.tissues';
fva_results.active = fva_Active;
fva_results.max = fva_Max;
fva_results.min = fva_Min;
save('checkpoints/fva_results.mat', 'fva_results')

%% test celfie tasks on constrained models
taskReportsCellfie = cell(numSamp, 1);
taskOkCellfie = zeros(length(tasks_Cellfie), numSamp);
for i = 1:numSamp 
    disp(['Model: ' num2str(i) ' of ' num2str(numSamp)])
        taskReportsCellfie{i} = checkTasks(closeModel(cfm_models_constrained{i}),[],false,false,false,tasks_Cellfie); 
        taskOkCellfie(:,i) = taskReportsCellfie{i}.ok;
end
saveTable = array2table(taskOkCellfie, 'VariableNames', expr_struct.tissues);
saveTable = [table({tasks_Cellfie.description}.', 'VariableNames', {'Task'}), saveTable];
writetable(saveTable, sprintf('%s/metabolic-tasks/cellfie_constrained.csv', directory_output), 'WriteRowNames', false);

%% save composition information for both constrained and unconstrained cfm models
model = ihuman; 
composition_cfm_models.rxns  = zeros([length(model.rxns), numSamp]);
composition_cfm_models.mets  = zeros([length(model.mets), numSamp]); 
composition_cfm_models_constrained.rxns  = zeros([length(model.rxns), numSamp]);
composition_cfm_models_constrained.mets  = zeros([length(model.mets), numSamp]); 
for i = 1:numSamp 
    % unconstrained models
    [~, irxn, ~] = intersect(model.rxns, cfm_models{i}.rxns, 'stable');
    [~, imet, ~] = intersect(model.mets, cfm_models{i}.mets, 'stable');
    composition_cfm_models.rxns(irxn, i) = 1;
    composition_cfm_models.mets(imet, i) = 1;
    % constrained
    [~, irxn, ~] = intersect(model.rxns, cfm_models_constrained{i}.rxns, 'stable');
    [~, imet, ~] = intersect(model.mets, cfm_models_constrained{i}.mets, 'stable');
    composition_cfm_models_constrained.rxns(irxn, i) = 1;
    composition_cfm_models_constrained.mets(imet, i) = 1;
    clear irxn imet
end

% humangem info
infoTable = table();
infoTable.rxnid     = transpose(1:length(model.rxns));
infoTable.rxnNames  = model.rxnNames;
infoTable.subSystems = model.subSystems;
infoTable.grRules   = model.grRules;
infoTable.eccodes   = model.eccodes;
infoTable.Properties.RowNames = model.rxns;
writetable(infoTable, sprintf('%s/composition_info.csv', directory_output), 'WriteRowNames', true);

% reactions
rxnsTable = [];
rxnsTable = array2table(composition_cfm_models.rxns);
rxnsTable.Properties.RowNames = model.rxns;
rxnsTable.Properties.VariableNames = expr_struct.tissues;
writetable(rxnsTable, sprintf('%s/composition_rxns.csv', directory_output), 'WriteRowNames', true);
rxnsTable = [];
rxnsTable = array2table(composition_cfm_models_constrained.rxns);
rxnsTable.Properties.RowNames = model.rxns;
rxnsTable.Properties.VariableNames = expr_struct.tissues;
writetable(rxnsTable, sprintf('%s/composition_rxns_constrained.csv', directory_output), 'WriteRowNames', true);

% metabolites
metsTable = [];
metsTable = array2table(composition_cfm_models.mets);
metsTable.Properties.RowNames = model.mets;
metsTable.Properties.VariableNames = expr_struct.tissues;
writetable(metsTable, sprintf('%s/composition_mets.csv', directory_output), 'WriteRowNames', true);
metsTable = [];
metsTable = array2table(composition_cfm_models_constrained.mets);
metsTable.Properties.RowNames = model.mets;
metsTable.Properties.VariableNames = expr_struct.tissues;
writetable(metsTable, sprintf('%s/composition_mets_constrained.csv', directory_output), 'WriteRowNames', true);

clear infoTable rxnsTable metsTable i


%% run simulations on unconstrained models
cfm_models_cobra   = cell(numSamp, 1);
pfba_Models   = cell(numSamp, 1);
pfba_Solutions = cell(numSamp, 1); 
knockout_Genes = struct();
knockout_Genes.grRatio   = cell(numSamp, 1);
knockout_Genes.grRateKO    = cell(numSamp, 1);
knockout_Genes.grRateWT    = cell(numSamp, 1);
knockout_Genes.hasEffect    = cell(numSamp, 1);
knockout_Genes.delRxns    = cell(numSamp, 1);
knockout_Genes.fluxSolution    = cell(numSamp, 1);
knockout_Reactions = knockout_Genes;
for i = 1:numSamp  
    disp(['Model: ' num2str(i) ' of ' num2str(numSamp)]) 
    % switch from RAVEN to COBRA
    curModel = cfm_models{i};
    curModel = ravenCobraWrapper(curModel);
    cfm_models_cobra{i} = curModel;

    % run pFBA and return further constrained models
    [~, ~, pfba_Models{i}, pfba_Solutions{i}] = pFBA(curModel, 'skipclass', 1);
    % save pfba solution
    saveTable = [];
    saveTable = array2table(pfba_Solutions{i}.x, 'VariableNames', {expr_struct.tissues{i}});
    saveTable = [table(pfba_Models{i}.rxns, 'VariableNames', {'rxn'}), saveTable];
    writetable(saveTable, sprintf('%s/pfba/%s.csv', directory_output, expr_struct.tissues{i}), 'WriteRowNames', false);
 
    % perform essentiality analysis on genes and reactions
    [knockout_Genes.grRatio{i}, knockout_Genes.grRateKO{i}, knockout_Genes.grRateWT{i}, knockout_Genes.hasEffect{i}, knockout_Genes.delRxns{i}, knockout_Genes.fluxSolution{i}] = singleGeneDeletion(curModel);
    [knockout_Reactions.grRatio{i}, knockout_Reactions.grRateKO{i}, knockout_Reactions.grRateWT{i}, knockout_Reactions.hasEffect{i}, knockout_Reactions.delRxns{i}, knockout_Reactions.fluxSolution{i}] = singleRxnDeletion(curModel);
    % save knockout results
    saveTable = [];
    saveTable = array2table(knockout_Reactions.grRatio{i}, 'VariableNames', {expr_struct.tissues{i}});
    saveTable = [table(knockout_Reactions.delRxns{i}, 'VariableNames', {'deleted rxn'}), saveTable];
    writetable(saveTable, sprintf('%s/knockout-reactions/%s.csv', directory_output, expr_struct.tissues{i}), 'WriteRowNames', false);
    saveTable = [];
    saveTable = array2table(knockout_Genes.grRatio{i}, 'VariableNames', {expr_struct.tissues{i}});
    saveTable = [table(curModel.geneNames, 'VariableNames', {'deleted gene'}), saveTable];
    writetable(saveTable, sprintf('%s/knockout-genes/%s.csv', directory_output, expr_struct.tissues{i}), 'WriteRowNames', false);

end
% save cobra models to save time for future analyses
save('checkpoints/cfm_models_cobra.mat', 'cfm_models_cobra')

% save pfba results
pfba_solutions.names = expr_struct.tissues';
pfba_solutions.models = pfba_Models;
pfba_solutions.solution = pfba_Solutions;
save('checkpoints/pfba_solutions.mat', 'pfba_solutions')

% save gene knockout results
knockout_simulations.names = expr_struct.tissues';
knockout_simulations.genes = knockout_Genes;
knockout_simulations.reactions = knockout_Reactions;
save('checkpoints/knockout_simulations.mat', 'knockout_simulations', '-v7.3')


%% run simulations on constrained models
cfm_models_constrained_cobra   = cell(numSamp, 1);
pfba_Models   = cell(numSamp, 1);
pfba_Solutions = cell(numSamp, 1); 
knockout_Genes = struct();
knockout_Genes.grRatio   = cell(numSamp, 1);
knockout_Genes.grRateKO    = cell(numSamp, 1);
knockout_Genes.grRateWT    = cell(numSamp, 1);
knockout_Genes.hasEffect    = cell(numSamp, 1);
knockout_Genes.delRxns    = cell(numSamp, 1);
knockout_Genes.fluxSolution    = cell(numSamp, 1);
knockout_Reactions = knockout_Genes;
for i = 1:numSamp  
    disp(['Model: ' num2str(i) ' of ' num2str(numSamp)]) 
    % switch from RAVEN to COBRA
    curModel = cfm_models_constrained{i};
    curModel = ravenCobraWrapper(curModel);
    cfm_models_constrained_cobra{i} = curModel;

    % run pFBA and return further constrained models
    [~, ~, pfba_Models{i}, pfba_Solutions{i}] = pFBA(curModel, 'skipclass', 1);
    % save pfba solution
    saveTable = [];
    saveTable = array2table(pfba_Solutions{i}.x, 'VariableNames', {expr_struct.tissues{i}});
    saveTable = [table(pfba_Models{i}.rxns, 'VariableNames', {'rxn'}), saveTable];
    writetable(saveTable, sprintf('%s/pfba-constrained/%s.csv', directory_output, expr_struct.tissues{i}), 'WriteRowNames', false);
 

    % perform essentiality analysis on genes and reactions
    [knockout_Genes.grRatio{i}, knockout_Genes.grRateKO{i}, knockout_Genes.grRateWT{i}, knockout_Genes.hasEffect{i}, knockout_Genes.delRxns{i}, knockout_Genes.fluxSolution{i}] = singleGeneDeletion(curModel);
    [knockout_Reactions.grRatio{i}, knockout_Reactions.grRateKO{i}, knockout_Reactions.grRateWT{i}, knockout_Reactions.hasEffect{i}, knockout_Reactions.delRxns{i}, knockout_Reactions.fluxSolution{i}] = singleRxnDeletion(curModel);
    % save knockout results
    saveTable = [];
    saveTable = array2table(knockout_Reactions.grRatio{i}, 'VariableNames', {expr_struct.tissues{i}});
    saveTable = [table(knockout_Reactions.delRxns{i}, 'VariableNames', {'deleted rxn'}), saveTable];
    writetable(saveTable, sprintf('%s/knockout-reactions-constrained/%s.csv', directory_output, expr_struct.tissues{i}), 'WriteRowNames', false);
    saveTable = [];
    saveTable = array2table(knockout_Genes.grRatio{i}, 'VariableNames', {expr_struct.tissues{i}});
    saveTable = [table(curModel.geneNames, 'VariableNames', {'deleted gene'}), saveTable];
    writetable(saveTable, sprintf('%s/knockout-genes-constrained/%s.csv', directory_output, expr_struct.tissues{i}), 'WriteRowNames', false);

end
% save cobra models to save time for future analyses
save('checkpoints/cfm_models_constrained_cobra.mat', 'cfm_models_constrained_cobra')

% save pfba results
pfba_solutions_constrained.names = expr_struct.tissues';
pfba_solutions_constrained.models = pfba_Models;
pfba_solutions_constrained.solution = pfba_Solutions;
save('checkpoints/pfba_solutions_constrained.mat', 'pfba_solutions_constrained')

% save gene knockout results
knockout_simulations_constrained.names = expr_struct.tissues';
knockout_simulations_constrained.genes = knockout_Genes;
knockout_simulations_constrained.reactions = knockout_Reactions;
save('checkpoints/knockout_simulations_constrained.mat', 'knockout_simulations_constrained', '-v7.3')


%% apply enzymatic constraints individually to each cfm model
% use default values
ModelAdapterManager.setDefault(fullfile(findGECKOroot,'tutorials','light_ecModel','HumanGEMAdapter.m')); 

% run GECKO pipeline
cfm_models_constrained_ec = cell(numSamp, 1);
for i = 1:numSamp  
    disp(['Model: ' num2str(i) ' of ' num2str(numSamp)])
    
    curModel = cfm_models_constrained{i};

    curecModel = makeEcModel(curModel, false); 
    curecModel = getECfromGEM(curecModel);
    kcatList   = fuzzyKcatMatching(curecModel);
    curecModel = selectKcatValue(curecModel, kcatList);
    curecModel = applyKcatConstraints(curecModel);
    curecModel = setProtPoolSize(curecModel);

    cfm_models_constrained_ec{i} = curecModel;
    cfm_models_constrained_ec{i}.id = curModel.id;

end
save('checkpoints/cfm_models_constrained_ec.mat', 'cfm_models_constrained_ec', '-v7.3')


%% run simulations on ec models
cfm_models_cobra_ec   = cell(numSamp, 1);
pfba_Models   = cell(numSamp, 1);
pfba_Solutions = cell(numSamp, 1); 
knockout_Genes = struct();
knockout_Genes.grRatio   = cell(numSamp, 1);
knockout_Genes.grRateKO    = cell(numSamp, 1);
knockout_Genes.grRateWT    = cell(numSamp, 1);
knockout_Genes.hasEffect    = cell(numSamp, 1);
knockout_Genes.delRxns    = cell(numSamp, 1);
knockout_Genes.fluxSolution    = cell(numSamp, 1);
knockout_Reactions = knockout_Genes;
for i = 1:numSamp  
    disp(['Model: ' num2str(i) ' of ' num2str(numSamp)]) 
    % switch from RAVEN to COBRA
    curModel = cfm_models_constrained_ec{i};
    curModel = ravenCobraWrapper(curModel);
    cfm_models_cobra_ec{i} = curModel;

    % run pFBA and return further constrained models
    [~, ~, pfba_Models{i}, pfba_Solutions{i}] = pFBA(curModel, 'skipclass', 1);
    % save pfba solution
    saveTable = [];
    saveTable = array2table(pfba_Solutions{i}.x, 'VariableNames', {expr_struct.tissues{i}});
    saveTable = [table(pfba_Models{i}.rxns, 'VariableNames', {'rxn'}), saveTable];
    writetable(saveTable, sprintf('%s/pfba-constrained-ec/%s.csv', directory_output, expr_struct.tissues{i}), 'WriteRowNames', false);
 
    % perform essentiality analysis on genes and reactions
    [knockout_Genes.grRatio{i}, knockout_Genes.grRateKO{i}, knockout_Genes.grRateWT{i}, knockout_Genes.hasEffect{i}, knockout_Genes.delRxns{i}, knockout_Genes.fluxSolution{i}] = singleGeneDeletion(curModel);
    %% THESE ARE NOT SAVED

    [knockout_Reactions.grRatio{i}, knockout_Reactions.grRateKO{i}, knockout_Reactions.grRateWT{i}, knockout_Reactions.hasEffect{i}, knockout_Reactions.delRxns{i}, knockout_Reactions.fluxSolution{i}] = singleRxnDeletion(curModel);
    % save knockout results
    saveTable = [];
    saveTable = array2table(knockout_Reactions.grRatio{i}, 'VariableNames', {expr_struct.tissues{i}});
    saveTable = [table(knockout_Reactions.delRxns{i}, 'VariableNames', {'deleted rxn'}), saveTable];
    writetable(saveTable, sprintf('%s/knockout-reactions-constrained-ec/%s.csv', directory_output, expr_struct.tissues{i}), 'WriteRowNames', false);
end

% save cobra models to save time for future analyses
save('checkpoints/cfm_models_cobra_ec.mat', 'cfm_models_cobra_ec')

% save pfba results
pfba_solutions_constrained_ec.names = expr_struct.tissues';
pfba_solutions_constrained_ec.models = pfba_Models;
pfba_solutions_constrained_ec.solution = pfba_Solutions;
save('checkpoints/pfba_solutions_constrained_ec.mat', 'pfba_solutions_constrained_ec')

% save gene knockout results
knockout_simulations_constrained_ec.names = expr_struct.tissues';
knockout_simulations_constrained_ec.genes = knockout_Genes;
knockout_simulations_constrained_ec.reactions = knockout_Reactions;
save('checkpoints/knockout_simulations_constrained_ec.mat', 'knockout_simulations_constrained_ec', '-v7.3')

%% Find rate limiting enzymes
rateLimitingEnzymes = cell(numSamp, 1);
rateLimitingEnzymesTable = table('Size', [0 4],'VariableTypes', {'string','string','string','string'}, 'VariableNames', {'enzyme', 'rxns', 'rxnNames', 'model'});
for i = 1:numSamp  
    disp(['Model: ' num2str(i) ' of ' num2str(numSamp)])
    
    curModel = cfm_models_constrained_ec{i};
    [~, rateLimitingEnzymes{i}] = sensitivityTuning(curModel);
    curTable = {};
    curTable(:,1) = rateLimitingEnzymes{i}.enzymes;
    curTable(:,2) = rateLimitingEnzymes{i}.rxns;
    curTable(:,3) = rateLimitingEnzymes{i}.rxnNames;
    curTable(:,4) = repmat({curModel.id}, length(curTable(:,3)), 1);
    rateLimitingEnzymesTable = [rateLimitingEnzymesTable; curTable];

end
writetable(rateLimitingEnzymesTable, sprintf('%s/other-simulations/ecmodels_rate_limiting_enzymes.csv', directory_output), 'WriteRowNames', false);


%% Sample models using RHMC
options.nStepsPerPoint = 100;   % sampling density
options.nPointsReturned = 1000; % number of points returned
for i = 6:numSamp  
    disp(['Model: ' num2str(i) ' of ' num2str(numSamp)])
    curModel = cfm_models_constrained{i}; 

    % sampling is time consuming, skip these.
    if curModel.id == "MARCO+ Macrophage (N)" || curModel.id == "MFAP5+ myofibroblasts (T)"
        disp('skipped.')
        continue
    end

    [~, samples_rhmc] =  sampleCbModel(curModel, [], 'RHMC', options);
    saveTable = array2table(samples_rhmc);
    saveTable.Properties.RowNames = curModel.rxns;
    writetable(saveTable, sprintf('%s/sampling-constrained/%s.csv', directory_output, expr_struct.tissues{i}), 'WriteRowNames', true);
end
