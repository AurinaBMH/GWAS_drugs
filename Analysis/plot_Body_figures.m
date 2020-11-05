% supplementary
clear all; close all; 
params = SetDefaultParams();
similarityTypes = {'MAGMAdefault', 'PPI_mapped_th600', 'eQTLbrain', 'AlleneQTLbrain'};
similarityTypes_label = {'SNP position', 'PPI network', 'Brain eQTL', 'AHBA'}; 

whatDiseases_GWAS = {'IBD','RA', 'HF', 'DIABETES'}; %'BIPandSCZ'
numGWAS = length(whatDiseases_GWAS); 
whatMeasures = 'allBody'; 
numDrugs = length(params.whatDiseases_Treatment); 
whatNull = sprintf('randomDrugR_%s_drugbank', params.whatTargets); 

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
    
    [rhosALL ,pValsALL] = DistinguishingCharBar(similarityTypes{s},whatProperty, whatNull, 'BF', whatDiseases_GWAS, true, length(similarityTypes), whatMeasures);
    figureName = sprintf('figures/BarChart_%s_%s', similarityTypes{s}, whatMeasures);
    print(gcf,figureName,'-dpng','-r300');
    
    
end




% does combinig scores improve matches?
% does combinig scores improve matches?
DOrecalc = true; 
f = plot_compareMeasures(whatDiseases_GWAS, whatMeasures, DOrecalc); 



