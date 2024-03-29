% Generating results for psychiatric disorders:
function [T_diabetes, T_bip] = plot_Psych_figures()
%-------------------------------------------------------
% Set the options
%-------------------------------------------------------
params = SetDefaultParams();
similarityTypes = {'MAGMAdefault', 'PPI_mapped_th600', 'eQTLbrain', 'AlleneQTLbrain'};
similarityTypes_label = {'SNP position', 'PPI network', 'Brain eQTL', 'AHBA'}; 
whatDiseases_GWAS = {'ADHD3','MDD4','SCZ3','BIP3','DIABETES2'};
numDrugs = length(params.whatDiseases_Treatment); 
whatMeasures = 'allPsych';
whatNull = sprintf('randomDrugR_%s_drugbank', params.whatTargets); 
numGWAS = length(whatDiseases_GWAS); 
V = nan(length(similarityTypes), numGWAS); 
whatDiseases_GWAS_name = cell(length(whatDiseases_GWAS), 1); 

% give names without numbers
for i=1:length(whatDiseases_GWAS) 
    name = whatDiseases_GWAS{i};
    whatDiseases_GWAS_name{i} = name(isstrprop(name,'alpha')); 
end

%-------------------------------------------------------
% Figure S1: Calculate pairwise matches 
%-------------------------------------------------------

nullScores = cell(length(similarityTypes),1); 

for s=1:length(similarityTypes)
    
    if contains(similarityTypes{s},'PPI')
        whatProperty = 'percPPIneighbors1';
    else
        if ~contains(similarityTypes{s},'Allen')
            whatProperty = 'P';
        elseif contains(similarityTypes{s},'Allen')
            whatProperty = 'zval';
        end
    end
    
    [rhosALL ,pValsALL, whatDiseases_Treatment, ~, enrichment_score_GWAS, enrichment_score_drug] = DistinguishingCharBar(similarityTypes{s},whatProperty, whatNull, 'BF', whatDiseases_GWAS, true, numDrugs, whatMeasures);
    % save scores for enrichment
    save(sprintf('enrichment_2024/enrichment_GWAS_%s.mat', similarityTypes{s}), 'enrichment_score_GWAS'); 
    save(sprintf('enrichment_2024/enrichment_drug_%s.mat', similarityTypes{s}), 'enrichment_score_drug')
    
    
    % find corresponsing match
    [T, INDr, INDc] = intersect(whatDiseases_Treatment, whatDiseases_GWAS_name, 'stable'); 
    % select disorder to itself - diagonal
    Pmatrix(s,:) = diag(pValsALL(INDr, INDc)); 

    figureName = sprintf('figures_2024/BarChart_psych_%s_%s_%s', similarityTypes{s}, whatMeasures, whatNull);
    print(gcf,figureName,'-dpng','-r300');

    
end

%-------------------------------------------------------
% Figure 2: Calculate matches for individual disorders
%-------------------------------------------------------
f = plot_measureOverview(Pmatrix, T, similarityTypes_label); 
figureName = sprintf('figures_2024/BarP_withinDisorder_%s_%s', whatMeasures, whatNull);
print(gcf,figureName,'-dpng','-r300');

% score genes by contribution: 
[Prank_diabetes, Drank_diabetes, Grank_diabetes] = rank_gene_contribution('DIABETES2', 'DIABETES', 'PPI_mapped_th600');
[Prank_bip, Drank_bip, Grank_bip] = rank_gene_contribution('BIP3', 'BIP', 'PPI_mapped_th600');

% these are genes that more significantly contribute to the match between
% GWAS and drugs compared to random; Prank_diabetes, Prank_bip; 
T_diabetes = join(Prank_diabetes, Drank_diabetes);
T_diabetes = sortrows(T_diabetes, 2, 'ascend');

T_bip = join(Prank_bip, Drank_bip);
T_bip = sortrows(T_bip, 2, 'ascend');

%-------------------------------------------------------
% Figure 3: Correspondence across different data processing methods
%-------------------------------------------------------
DOrecalc = true; % select true if any data was updated since the last run
pPlot_all = plot_compareMeasures(whatDiseases_GWAS, whatMeasures, DOrecalc); 

%-------------------------------------------------------
% Figure S2: Correlation between GWAS-based scores for different data processing methods
%-------------------------------------------------------
Mnames = cell(numGWAS,1); 
Mnumbers = cell(numGWAS,1); 

for i=1:numGWAS

    [f, Mnames{i}, Mnumbers{i}] = correlate_geneMeasures(whatDiseases_GWAS{i}, whatMeasures, true);
    figureName = sprintf('figures_2024/%s_geneMeasures_%s', whatDiseases_GWAS{i}, whatMeasures);
    print(f,figureName,'-dpng','-r300');
    
end

%-------------------------------------------------------
% Figure S4: different null distributions for bipolar disorder and diabetes
%-------------------------------------------------------
% Plot null distributions when choosing from all and from psychiatric drugs: 
% in this example PPI-based significant results are used: bipolar disorder and diabetes; 
[f, Ptable_null] = plot_nullDistributions();
figureName = 'figures_2024/Null_distribution_comparison';
print(f,figureName,'-dpng','-r300');

end


    
    