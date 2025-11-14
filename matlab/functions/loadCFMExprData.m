function [expr_data] = loadCFMExprData(filepath)
%  Loads and processes patient metadata and gene expression data
%  and returns a structure containing processed metadata and gene expression levels.
%
%  Input:
%      filename - String, name of the CSV file containing patient metadata
%
%  Output:
%      expr_data - Structure containing processed data fields


% filename = 'ellison-patients-metadata.csv'

% load expression profiles (obtained with r)
expressions = readtable(filepath);

% fix column names
name_mapping = {   'Var1', 'Genes';
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
expr_data.genes = expressions.Genes; % gene names in ENSG format
expr_data.clusters = expressions.Properties.VariableNames(2:end)'; % celltype names

expr_data.maintypes = [repmat("Cancer", 1, 2),  ...
                                    repmat("Fibroblast", 1, 14), ...
                                    repmat("Macrophage", 1, 6)]'; % main type

expr_data.tissues = cell(size(expr_data.clusters));
for i = 1:length(expr_data.clusters)
    if expr_data.clusters{i}(end-1) == 'T'
        expr_data.tissues{i} = 'Tumor';
    elseif expr_data.clusters{i}(end-1) == 'N'
        expr_data.tissues{i} = 'Normal';
    end
end

expr_data.cpm = table2array(expressions(:, 2:end)); % cpm values 
expr_data.cpm_scaled = expr_data.cpm .* (1e6 ./ sum(expr_data.cpm)); % rescaled cpm


