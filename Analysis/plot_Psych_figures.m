% matching for psychiatric disorders
clear all; close all; 

params = SetDefaultParams();
similarityTypes = {'MAGMAdefault', 'PPI_mapped_th600', 'eQTLbrain', 'AlleneQTLbrain'};
similarityTypes_label = {'SNP position', 'PPI network', 'Brain eQTL', 'AHBA'}; 
whatDiseases_GWAS = {'ADHD2', 'MDD3','SCZ3','BIP3','DIABETES'};
numDrugs = length(params.whatDiseases_Treatment); 
whatMeasures = 'allPsych';
whatNull = sprintf('randomDrugR_%s_drugbank', params.whatTargets); 


numGWAS = length(whatDiseases_GWAS); 
V = nan(length(similarityTypes), numGWAS); 

% give names without numbers
for i=1:length(whatDiseases_GWAS) 
    name = whatDiseases_GWAS{i};
    whatDiseases_GWAS_name{i} = name(isstrprop(name,'alpha')); 
end

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
    save(sprintf('enrichment_2022/enrichment_GWAS_%s.mat', similarityTypes{s}), 'enrichment_score_GWAS'); 
    save(sprintf('enrichment_2022/enrichment_drug_%s.mat', similarityTypes{s}), 'enrichment_score_drug')
    
    
    % find corresponsing match
    [T_diabetes, INDr, INDc] = intersect(whatDiseases_Treatment, whatDiseases_GWAS_name, 'stable'); 
    % select disorder to itself - diagonal
    Pmatrix(s,:) = diag(pValsALL(INDr, INDc)); 

    figureName = sprintf('figures_2022/BarChart_%s_%s_%s', similarityTypes{s}, whatMeasures, whatNull);
    print(gcf,figureName,'-dpng','-r300');
    % save scores for enrichment
    
    
end

f = plot_measureOverview(Pmatrix, T_diabetes, similarityTypes_label); 
figureName = sprintf('figures_2022/BarP_withinDisorder_%s_%s', whatMeasures, whatNull);
print(gcf,figureName,'-dpng','-r300');

% score genes by contribution: 
[Prank_diabetes, Drank_diabetes, Grank_diabetes] = rank_gene_contribution('DIABETES', 'DIABETES', 'PPI_mapped_th600');
[Prank_bip, Drank_bip, Grank_bip] = rank_gene_contribution('BIP2', 'BIP', 'PPI_mapped_th600');

T_diabetes = join(Drank_diabetes,Grank_diabetes);
T_diabetes = sortrows(T_diabetes, 3, 'descend');

T_bip = join(Drank_bip,Grank_bip);
T_bip = sortrows(T_bip, 3, 'descend');

% plot null distributions when choosing from all and from psychiatric
% drugs: in this example: use PPI-significant results: BIP GWAS vs BIP drugs and DIABETES vs DIABETES drugs
[f, Ptable_null] = plot_nullDistributions(); 
figureName = 'figures_2022/Null_distribution_comparison';
print(f,figureName,'-dpng','-r300');


% plot correlation between different measures for one representative disorder
for i=1:numGWAS

    [f, Mnames{i}, Mnumbers{i}] = correlate_geneMeasures(whatDiseases_GWAS{i}, whatMeasures, true);
    figureName = sprintf('figures_2022/%s_geneMeasures_%s', whatDiseases_GWAS{i}, whatMeasures);
    print(f,figureName,'-dpng','-r300');
    
end

% does combinig scores improve matches?
DOrecalc = true; 
f = plot_compareMeasures(whatDiseases_GWAS, whatMeasures, DOrecalc); 








%% don't plot these - will say in words, it's a null result
% f = figure('color','w', 'Position', [300, 300, 1500, 400]); 
% % choose more diverging colours
% barColor = BF_getcmap('set5',5);  
% 
% for i=1:numGWAS
%     
%     ax{i} = subplot(1,numGWAS,i); hold on
%     title(sprintf('%s', whatDiseases_GWAS{i}))
%     b = bar(Ptable.(whatDiseases_GWAS{i}).Pvals);
%     ylabel('-log10(P)')
%     set(gca,'FontSize', 14)
%     b.CData = barColor;
%     b.FaceColor = 'flat';
%     
%     ax{i}.XTick = 1:length(similarityTypes)+1;
%     ax{i}.XTickLabel = [similarityTypes_label, 'Combined'];
%     ax{i}.XTickLabelRotation = 45;
%     
% end
% 
% % rescale axes
% linkaxes([ax{:}],'y');
% 
% print(gcf,'figures/comparePvals','-dpng','-r300');

    
    