function [f, nr_g] = count_number_contributing_genes()

% count the number of contributing genes for each disorder and mapping
% method for revision. 
params = SetDefaultParams(); 
[geneNamesDrug,drugScoresAll] = GiveMeNormalizedScoreVectors(params.whatDiseases_Treatment_ALL,'Drug');

disorders = {'ADHD3','BIP3','SCZ3', 'MDD4','DIABETES2'}; 
nr_g = zeros(length(disorders),4); 

for i=1:length(disorders)
       
ix_drug = find(drugScoresAll(:,i));
load(sprintf('resultsTable_%s_BF_2024_all_drugbank.mat', disorders{i}))


SCZ_gwas = geneScores.MAGMAdefault.ZSTAT; 
SCZ_gwas(isnan(SCZ_gwas)) = 0; 
ix_gwas = find(SCZ_gwas); 
nr_g(i,1) = length(intersect(ix_gwas, ix_drug));
fprintf('MAGMA %d genes contribute to similarity score\n', nr_g(i,1))

SCZ_gwas = geneScores.PPI_mapped_th600.numPPIneighbors1; 
SCZ_gwas(isnan(SCZ_gwas)) = 0; 
ix_gwas = find(SCZ_gwas); 
nr_g(i,2) = length(intersect(ix_gwas, ix_drug));
fprintf('PPI %d genes contribute to similarity score\n', nr_g(i,2))

SCZ_gwas = geneScores.eQTLbrain.ZSTAT; 
SCZ_gwas(isnan(SCZ_gwas)) = 0; 
ix_gwas = find(SCZ_gwas); 
nr_g(i,3) = length(intersect(ix_gwas, ix_drug));
fprintf('eQTL %d genes contribute to similarity score\n', nr_g(i,3))

SCZ_gwas = geneScores.AllenMapped.zval; 
SCZ_gwas(isnan(SCZ_gwas)) = 0; 
ix_gwas = find(SCZ_gwas); 
nr_g(i,4) = length(intersect(ix_gwas, ix_drug));
fprintf('AHBA %d genes contribute to similarity score\n', nr_g(i,4))

end

f = figure('color','w', 'Position', [300, 300, 800, 800]);

% plot
xvalues = {'MAGMAdefault', 'PPI mapped th600', 'eQTLbrain', 'AHBA'}; 
yvalues = {'ADHD', 'BIP', 'SCZ', 'MDD', 'DIABETES'};
h = heatmap(xvalues,yvalues,nr_g);
colors = cbrewer('seq', 'Reds', 64);
colormap(colors); 
set(gca,'FontSize', 18)
colorbar
title('Number of contributing genes')

end









