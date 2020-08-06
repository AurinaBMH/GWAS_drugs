clear all; 

numNulls = 1000; 
whatDiseases_Treatment = {'ADHD','BIP','SCZ','MDD','pulmonary','cardiology','gastro','diabetes'};
params = SetDefaultParams();

[geneNamesGWAS,geneWeightsGWAS] = GiveMeNormalizedScoreVector('ADHD','GWAS','MAGMAdefault',params.geneScore, params.whatThreshold);
[~,~, disorderDrugs, allDrugs] = ImportTreatmentLists(true, params.whatDrugTargets);

numGenes = length(geneNamesGWAS); 

RNADOMdrugs = zeros(numGenes, numNulls); 
RNADOMdrugs_treatment = cell(length(whatDiseases_Treatment),1); 

for i = 1:length(whatDiseases_Treatment)
    for j = 1:numNulls
        diseaseName = whatDiseases_Treatment{i}; 
        RNADOMdrugs(:,j) = give_randomDrug_null(diseaseName, disorderDrugs, allDrugs);
    end
    RNADOMdrugs_treatment{i} = RNADOMdrugs;
end

fileName = sprintf('DataOutput/nulls_%d_randomDrug.mat', numNulls); 
save(fileName, 'RNADOMdrugs_treatment', 'whatDiseases_Treatment', 'params'); 