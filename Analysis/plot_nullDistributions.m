function plot_nullDistributions()

whatMeasures = 'allPsych';
colPsy = [253,141,6]/255; 
colAll = [115,115,115]/255; 
colLine = [252,78,42]/255; 

% for several borderline matches plot null distributions and real data
% based on all radnom drug nulls or only using psychiatric drug nulls; 
typesNull = {'randomDrugP','randomDrugP_psych'}; 
for t=1:length(typesNull)
% for SCZ
whatDiseases_GWAS = {'SCZ'}; 
[SCZ_rhosALL{t} ,~, whatDiseases_Treatment, ~, ~, SCZ_null{t}] = ...
    DistinguishingCharBar('PPI_mapped_th600','percPPIneighbors1', typesNull{t}, 'BF', whatDiseases_GWAS, false, 4, whatMeasures);

% for DIABETES
whatDiseases_GWAS = {'DIABETES'}; 
[DIABETES_rhosALL{t} ,~, whatDiseases_Treatment, ~, ~, DIABETES_null{t}] = ...
    DistinguishingCharBar('PPI_mapped_th600','percPPIneighbors1', typesNull{t}, 'BF', whatDiseases_GWAS, false, 4, whatMeasures);
end

f = figure('color','w', 'Position', [300, 300, 1200, 500]);

% find BIP in SCZ list
[~, IND_scz] = intersect(whatDiseases_Treatment,'BIP', 'stable');
rho_SCZ_BIP = SCZ_rhosALL{1}(IND_scz); 
N_scz_all = SCZ_null{1}{1}{IND_scz}; 
N_scz_psy = SCZ_null{2}{1}{IND_scz};

ax{1} = subplot(1,2,1); 
histogram(N_scz_all, 50, 'EdgeColor', colAll, 'FaceColor', [1 1 1], 'LineWidth', 2); hold on; 
histogram(N_scz_psy, 50, 'EdgeColor', colPsy, 'FaceColor', [1 1 1], 'LineWidth', 2); hold on; 
plot([rho_SCZ_BIP rho_SCZ_BIP],[0 300], 'LineWidth', 3, 'Color', colLine)
box off

xlabel('GWAS-treatment similarity')
ylabel({'Count'})
set(gca,'FontSize', 20)
title('Psychiatric disorder')
        
        

% find DIABETES in DIABETES list
[~, IND_diabetes] = intersect(whatDiseases_Treatment,'DIABETES', 'stable');
rho_DIABETES = DIABETES_rhosALL{1}(IND_diabetes); 
N_diabetes_all = DIABETES_null{1}{1}{IND_diabetes}; 
N_diabetes_psy = DIABETES_null{2}{1}{IND_diabetes};

ax{2} = subplot(1,2,2);
histogram(N_diabetes_all, 50, 'EdgeColor', colAll, 'FaceColor', [1 1 1], 'LineWidth', 2); hold on; 
histogram(N_diabetes_psy, 50, 'EdgeColor', colPsy, 'FaceColor', [1 1 1], 'LineWidth', 2); hold on; 
plot([rho_DIABETES rho_DIABETES],[0 300], 'LineWidth', 3, 'Color', colLine)
box off

xlabel('GWAS-treatment similarity')
ylabel({'Count'})
set(gca,'FontSize', 20)
title('Diabetes')

linkaxes([ax{:}],'y');
linkaxes([ax{:}],'x');

end