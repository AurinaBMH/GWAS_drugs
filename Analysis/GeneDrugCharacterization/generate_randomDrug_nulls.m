numNulls = 1000; 
whatDiseases_Treatment = {'ADHD','BIP','SCZ','MDD','pulmonary','cardiology','gastro','diabetes'};

RNADOMdrugs = zeros(numGenes, numNulls); 
RNADOMdrugs_treatment = cell(length(whatDiseases_Treatment),1); 

for i = 1:length(whatDiseases_Treatment)
    for j = 1:numNulls
        RNADOMdrugs(:,j) = give_randomDrug_null(diseaseName, disorderDrugs, allDrugs);
    end
    RNADOMdrugs_treatment{i} = RNADOMdrugs;
end