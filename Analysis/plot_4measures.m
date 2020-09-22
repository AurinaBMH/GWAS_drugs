% matching for psychiatric disorders
similarityTypes = {'MAGMAdefault', 'PPI_mapped_th900', 'eQTLbrain', 'AllenMeanCoexpMapped'}; 
whatDiseases_GWAS = {'ADHD','MDD2','SCZ','BIP2','DIABETES'}; 

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
    
    [rhosALL ,pValsALL] = DistinguishingCharBar(similarityTypes{s},whatProperty, 'randomDrugP', 'BF', whatDiseases_GWAS);
    
end