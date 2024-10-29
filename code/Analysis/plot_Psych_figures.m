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
% Figure S2: Calculate pairwise matches 
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
    
    [rhosALL ,pValsALL_randDrug, whatDiseases_Treatment_randDrug, null_scores_all, enrichment_score_GWAS_randDrug, enrichment_score_drug_randDrug] = ...
        DistinguishingCharBar(similarityTypes{s}, whatProperty, whatNull, 'BF', whatDiseases_GWAS, true, numDrugs, whatMeasures);
    % save scores for enrichment
    save(sprintf('enrichment_2024/enrichment_GWAS_%s.mat', similarityTypes{s}), 'enrichment_score_GWAS_randDrug'); 
    save(sprintf('enrichment_2024/enrichment_drug_%s.mat', similarityTypes{s}), 'enrichment_score_drug_randDrug')
    
    
    % find corresponsing match
    [T, INDr, INDc] = intersect(whatDiseases_Treatment_randDrug, whatDiseases_GWAS_name, 'stable'); 
    % select disorder to itself - diagonal
    Pmatrix(s,:) = diag(pValsALL_randDrug(INDr, INDc)); 
    rhomatrix = diag(rhosALL(INDr, INDc)); 

    % get p, vals, rho vals, mean rand rho and SD 
    P_all_order = pValsALL_randDrug(INDr, INDc); 
    RHO_all_order = rhosALL(INDr, INDc); 
    for ddd = 1:length(whatDiseases_GWAS)
        % reorder based on treatment order
        whatDisease = whatDiseases_GWAS{ddd}; 
        NULL = null_scores_all.(whatDisease); 
        NULL_ord = NULL(INDr); 
        NULL_all_order.(whatDisease) = NULL_ord; 
        
    end

    figureName = sprintf('figures_2024/BarChart_psych_%s_%s_%s', similarityTypes{s}, whatMeasures, whatNull);
    print(gcf,figureName,'-dpng','-r300');

    
    % get p-values for psych null
    [~ ,pValsALL_psychDrug, whatDiseases_Treatment_psychDrug] = ...
        DistinguishingCharBar(similarityTypes{s},whatProperty, 'randomDrugP_all_drugbank_psych', 'BF', whatDiseases_GWAS, false, numDrugs, whatMeasures);
     % find corresponsing match
    [T_psych, INDr_psych, INDc_psych] = intersect(whatDiseases_Treatment_psychDrug, whatDiseases_GWAS_name, 'stable'); 
    % select disorder to itself - diagonal
    Pmatrix_psych(s,:) = diag(pValsALL_psychDrug(INDr_psych, INDc_psych)); 
    
    
end

%-------------------------------------------------------
% Figure 2: Calculate matches for individual disorders
%-------------------------------------------------------
f = plot_measureOverview(Pmatrix, T, similarityTypes_label); 
figureName = sprintf('figures_2024/BarP_withinDisorder_%s_%s', whatMeasures, whatNull);
print(gcf,figureName,'-dpng','-r300');

% f = plot_measureOverview(Pmatrix_psych, T_psych, similarityTypes_label); 
% figureName = sprintf('figures_2024/BarP_withinDisorder_%s_%s', whatMeasures, 'randomDrugP_all_drugbank_psych');
% print(gcf,figureName,'-dpng','-r300');

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
DOrecalc = false; % select true if any data was updated since the last run
pPlot_all = plot_compareMeasures(whatDiseases_GWAS, whatMeasures, DOrecalc); 

%-------------------------------------------------------
% Figure S3: Correlation between GWAS-based scores for different data processing methods
%-------------------------------------------------------
Mnames = cell(numGWAS,1); 
Mnumbers = cell(numGWAS,1); 

for i=1:numGWAS

    [f, Mnames{i}, Mnumbers{i}] = correlate_geneMeasures(whatDiseases_GWAS{i}, whatMeasures, true);
    figureName = sprintf('figures_2024/%s_geneMeasures_%s', whatDiseases_GWAS{i}, whatMeasures);
    print(f,figureName,'-dpng','-r300');
    
end

%-------------------------------------------------------
% Figure S8: different null distributions for bipolar disorder and diabetes
%-------------------------------------------------------
% Plot null distributions when choosing from all and from psychiatric drugs: 
% in this example PPI-based significant results are used: bipolar disorder and diabetes; 
[f, Ptable_null] = plot_nullDistributions();
figureName = 'figures_2024/Null_distribution_comparison';
print(f,figureName,'-dpng','-r300');

%-------------------------------------------------------
% Figure in eMethods 4: The number of genes contributing to the similarity score for each disorder and mapping method.
%-------------------------------------------------------
% plot a matrix indicating a number of contributing genes in each disorder
% and mapping method
[f, nr_g] = count_number_contributing_genes(); 

figureName = 'figures_2024/Number_of_contributing_genes';
print(f,figureName,'-dpng','-r300');

end


    
    