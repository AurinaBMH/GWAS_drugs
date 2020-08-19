function generate_randomDrug_nulls()
 
% this function generates a set of random nulls by selecting the selected
% numnber of genes each time; Overall, the nulls for different numbers of
% drugs are very similar, could use only one null with the average number
% of drugs for comparison to all, but for now making a separate null for
% each drug list; 

numNulls = 5000; 
whatDiseases_Treatment = {'ADHD','BIP','SCZ','MDD','pulmonary','cardiology','gastro','diabetes'};
params = SetDefaultParams();

[geneNamesGWAS,~] = GiveMeNormalizedScoreVector('ADHD','GWAS','MAGMAdefault',params.geneScore, params.whatThreshold);
[~,~, disorderDrugs, allDrugs] = ImportTreatmentLists(true, params.whatDrugTargets);

numGenes = length(geneNamesGWAS); 

RANDOMdrugs = zeros(numGenes, numNulls); 
RANDOMdrugs_treatment = cell(length(whatDiseases_Treatment),1); 

for i = 1:length(whatDiseases_Treatment)
    for j = 1:numNulls
        diseaseName = whatDiseases_Treatment{i}; 
        RANDOMdrugs(:,j) = give_randomDrug_null(diseaseName, disorderDrugs, allDrugs);
    end
    RANDOMdrugs_treatment{i} = RANDOMdrugs;
end

fileName = sprintf('DataOutput/nulls_%d_%stargets_randomDrug.mat', numNulls, params.whatDrugTargets); 
save(fileName, 'RANDOMdrugs_treatment', 'whatDiseases_Treatment', 'params');
end