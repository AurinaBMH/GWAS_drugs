% reprocessing data using DrugBank data: 
dataTable = give_drugTargets('active', 'drugbank'); 
GenerateResultsTables

% generate nulls
generate_randomDrug_nulls('drugbank')

% generate results
% plot_Psych_figures.m
f = plot_measureOverview(Pmatrix, T, similarityTypes_label);


f = plot_nullDistributions();

correlate_geneMeasures(whatDiseases_GWAS{i}, whatMeasures, true);

f = plot_compareMeasures(whatDiseases_GWAS, whatMeasures, DOrecalc);