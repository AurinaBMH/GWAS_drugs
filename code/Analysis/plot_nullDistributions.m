function [f, Ptable] = plot_nullDistributions(whatYear)

if nargin < 1
    whatYear = '2024'; 
end

whatMeasures = 'allPsych';
colPsy = [253,141,6]/255; 
colAll = [115,115,115]/255; 
colLine = [252,78,42]/255; 

% for several borderline matches plot null distributions and real data
% based on all radnom drug nulls or only using psychiatric drug nulls; 
typesNull = {'randomDrugR_all_drugbank','randomDrugP_all_drugbank_psych'}; 
for t=1:length(typesNull)
    
% for BIP
if strcmp(whatYear, '2024')
    whatDiseases_GWAS = {'BIP3'};
elseif strcmp(whatYear, '2021')
    whatDiseases_GWAS = {'BIP2'};
end
[BIP_rhosALL{t} ,~, whatDiseases_Treatment, BIP_null{t}] = ...
    DistinguishingCharBar('PPI_mapped_th600','percPPIneighbors1', typesNull{t}, 'BF', whatDiseases_GWAS, false, 4, whatMeasures);

% for DIABETES, use newer GWAS
if strcmp(whatYear, '2024')
    whatDiseases_GWAS = {'DIABETES2'};
elseif strcmp(whatYear, '2021')
    whatDiseases_GWAS = {'DIABETES'};
end

[DIABETES_rhosALL{t} ,~, whatDiseases_Treatment, DIABETES_null{t}] = ...
    DistinguishingCharBar('PPI_mapped_th600','percPPIneighbors1', typesNull{t}, 'BF', whatDiseases_GWAS, false, 4, whatMeasures);
end

f = figure('color','w', 'Position', [300, 300, 1200, 500]);

% find BIP in BIP2 list
[~, IND_bip] = intersect(whatDiseases_Treatment,'BIP', 'stable');
rho_BIP_BIP = BIP_rhosALL{1}(IND_bip); 
N_bip_all = BIP_null{1}{1}{IND_bip}; 
N_bip_psy = BIP_null{2}{1}{IND_bip};

Ptable = struct; 
Ptable.BIPall = mean(N_bip_all>rho_BIP_BIP); 
Ptable.BIPpsy = mean(N_bip_psy>rho_BIP_BIP); 

ax{1} = subplot(1,2,1); 
histogram(N_bip_all, 50, 'EdgeColor', colAll, 'FaceColor', [1 1 1], 'LineWidth', 2); hold on; 
histogram(N_bip_psy, 50, 'EdgeColor', colPsy, 'FaceColor', [1 1 1], 'LineWidth', 2); hold on; 
plot([rho_BIP_BIP rho_BIP_BIP],[0 600], 'LineWidth', 3, 'Color', colLine)
box off
legend('All treatment null','Psychiatric treatment null')
legend boxoff  

xlabel('GWAS-treatment similarity')
ylabel({'Count'})
set(gca,'FontSize', 20)
title('Bipolar disorder')

% find DIABETES in DIABETES list
[~, IND_diabetes] = intersect(whatDiseases_Treatment,'DIABETES', 'stable');
rho_DIABETES = DIABETES_rhosALL{1}(IND_diabetes); 
N_diabetes_all = DIABETES_null{1}{1}{IND_diabetes}; 
N_diabetes_psy = DIABETES_null{2}{1}{IND_diabetes};

Ptable.DIABall = mean(N_diabetes_all>rho_DIABETES); 
Ptable.DIABpsy = mean(N_diabetes_psy>rho_DIABETES); 

ax{2} = subplot(1,2,2);
histogram(N_diabetes_all, 50, 'EdgeColor', colAll, 'FaceColor', [1 1 1], 'LineWidth', 2); hold on; 
histogram(N_diabetes_psy, 50, 'EdgeColor', colPsy, 'FaceColor', [1 1 1], 'LineWidth', 2); hold on; 
plot([rho_DIABETES rho_DIABETES],[0 600], 'LineWidth', 3, 'Color', colLine)
box off
legend('All treatment null','Psychiatric treatment null')
legend boxoff  

xlabel('GWAS-treatment similarity')
ylabel({'Count'})
set(gca,'FontSize', 20)
title('Diabetes')


linkaxes([ax{:}],'y');
linkaxes([ax{:}],'x');

end