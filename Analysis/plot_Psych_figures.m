% matching for psychiatric disorders
clear all; close all; 
params = SetDefaultParams();
similarityTypes = {'MAGMAdefault', 'PPI_mapped_th600', 'eQTLbrain', 'AllenMeanCoexpMapped'};
similarityTypes_label = {'SNP position', 'PPI network', 'Brain eQTL', 'AHBA'}; 
whatDiseases_GWAS = {'ADHD','MDD2','SCZ','BIP2','DIABETES'};
whatMeasures = 'allPsych';
whatNull = 'randomDrugP'; 


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
            whatProperty = 'r';
        end
    end
    
    [rhosALL ,pValsALL, whatDiseases_Treatment] = DistinguishingCharBar(similarityTypes{s},whatProperty, whatNull, 'BF', whatDiseases_GWAS, true, length(whatDiseases_GWAS)-1, whatMeasures);
    % find corresponsing match
    [T, INDr, INDc] = intersect(whatDiseases_Treatment, whatDiseases_GWAS_name, 'stable'); 
    % select disorder to itself - diagonal
    Pmatrix(s,:) = diag(pValsALL(INDr, INDc)); 

    figureName = sprintf('figures/BarChart_%s_%s', similarityTypes{s}, whatMeasures);
    print(gcf,figureName,'-dpng','-r300');
    
    
end

f = plot_measureOverview(Pmatrix, T, similarityTypes_label); 
figureName = sprintf('figures/BarP_withinDisorder_%s', whatMeasures);
print(gcf,figureName,'-dpng','-r300');


% for several borderline matches plot null distributions and real data
% based on all radnom drug nulls or only using psychiatric drug nulls; 
typesNull = {'randomDrugP','randomDrugP_psych'}; 
for t=1:length(typesNull)
% for SCZ
[SCZ_rhosALL{t} ,~, whatDiseases_Treatment, ~, ~, SCZ_null{t}] = ...
    DistinguishingCharBar('PPI_mapped_th600','percPPIneighbors1', typesNull{t}, 'BF', {'SCZ'}, true, length(whatDiseases_GWAS), whatMeasures);

% for DIABETES
[DIABETES_rhosALL{t} ,~, whatDiseases_Treatment, ~, ~, DIABETES_null{t}] = ...
    DistinguishingCharBar('PPI_mapped_th600','percPPIneighbors1', typesNull{t}, 'BF', {'DIABETES'}, true, length(whatDiseases_GWAS), whatMeasures);
end

% find BIP in SCZ list
[~, IND_scz] = intersect(whatDiseases_Treatment,'BIP', 'stable');
rho_SCZ_BIP = SCZ_rhosALL{1}(IND_scz); 
N_scz_all = SCZ_null{1}{1}{IND_scz}; 
N_scz_psy = SCZ_null{2}{1}{IND_scz};



% find DIABETES in DIABETES list
[~, IND_diabetes] = intersect(whatDiseases_Treatment,'DIABETES', 'stable');
rho_DIABETES = DIABETES_rhosALL{1}(IND_diabetes); 
N_diabetes_all = DIABETES_null{1}{1}{IND_diabetes}; 
N_diabetes_psy = DIABETES_null{2}{1}{IND_diabetes};




% plot correlation between different measures for one representative disorder
for i=1:numGWAS

    [f, Mnames{i}, Mnumbers{i}] = correlate_geneMeasures(whatDiseases_GWAS{i}, whatMeasures, true);
    figureName = sprintf('figures/%s_geneMeasures_%s', whatDiseases_GWAS{i}, whatMeasures);
    print(f,figureName,'-dpng','-r300');
    
end


% does combinig scores improve matches?
DOrecalc = false; 
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

    
    