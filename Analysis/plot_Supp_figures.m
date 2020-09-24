% supplementary
clear all; close all; 
similarityTypes = {'MAGMAdefault', 'PPI_mapped_th600', 'eQTLbrain', 'AllenMeanCoexpMapped'};
similarityTypes_label = {'SNP position', 'PPI network', 'Brain eQTL', 'AHBA'}; 

whatDiseases_GWAS = {'IBD','RA', 'HF', 'DIABETES'}; %'BIPandSCZ'
numGWAS = length(whatDiseases_GWAS); 
whatMeasures = 'allBody'; 

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
    
    [rhosALL ,pValsALL] = DistinguishingCharBar(similarityTypes{s},whatProperty, 'randomDrugP', 'BF', whatDiseases_GWAS, true, length(similarityTypes), whatMeasures);
    figureName = sprintf('figures/BarChart_%s_%s', similarityTypes{s}, whatMeasures);
    print(gcf,figureName,'-dpng','-r300');
    
    
end







% plot correlation between different measures for one representative disorder
%f = figure('color','w', 'Position', [300, 300, 2000, 1000]);
for i=1:numGWAS
    
    %ax{i} = subplot(1,numGWAS, i); hold on;
    [f, Mnames{i}, Mnumbers{i}] = correlate_geneMeasures(whatDiseases_GWAS{i}, whatMeasures, true);
    figureName = sprintf('figures/%s_geneMeasures_%s', whatDiseases_GWAS{i}, whatMeasures);
    print(f,figureName,'-dpng','-r300');
end


% does combinig scores improve matches?
plotHow = 'horizontal'; 
f = plot_compareMeasures(whatDiseases_GWAS, whatMeasures, true, plotHow); 
figureName = sprintf('figures/compareMeasures_%s_%s', whatMeasures, plotHow);
print(f,figureName,'-dpng','-r300');