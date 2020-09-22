% matching for psychiatric disorders
clear all; close all; 
similarityTypes = {'MAGMAdefault', 'PPI_mapped_th600', 'eQTLbrain', 'AllenMeanCoexpMapped'}; 
whatDiseases_GWAS = {'ADHD','MDD2','SCZ','BIP2','DIABETES'}; %'BIPandSCZ'

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
    
    [rhosALL ,pValsALL] = DistinguishingCharBar(similarityTypes{s},whatProperty, 'randomDrugP', 'BF', whatDiseases_GWAS, true, length(similarityTypes));
    figureName = sprintf('figures/BarChart_%s', similarityTypes{s});
    print(gcf,figureName,'-dpng','-r300');
end