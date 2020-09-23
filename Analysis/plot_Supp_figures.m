% supplementary
whatDiseases_GWAS = {'IBD','RA', 'HF', 'DIABETES'}; %'BIPandSCZ'
numGWAS = length(whatDiseases_GWAS); 
whatMeasures = 'allBody'; 
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