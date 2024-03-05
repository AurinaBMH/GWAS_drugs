function generate_randomDrug_nulls(whatSelection)
if nargin<1
    whatSelection = 'drugbank'; 
end
% this function generates a set of random nulls by selecting the selected
% numnber of genes each time; Overall, the nulls for different numbers of
% drugs are very similar, could use only one null with the average number
% of drugs for comparison to all, but for now making a separate null for
% each drug list; 

numNulls = 5000; 
params = SetDefaultParams();
whatDiseases_Treatment = params.whatDiseases_Treatment_ALL; 

[geneNamesGWAS,~] = GiveMeNormalizedScoreVector('ADHD3','GWAS','MAGMAdefault',params.geneScore, params.whatThreshold);
[~,~, disorderDrugs, allDrugs] = ImportTreatmentLists(true, params.whatDrugTargets, params.whatTargets);

numGenes = length(geneNamesGWAS); 
RANDOMdrugs = zeros(numGenes, numNulls); 
RANDOMdrugs_treatment = cell(length(whatDiseases_Treatment),1); 

for i = 1:length(whatDiseases_Treatment)
    for j = 1:numNulls
        diseaseName = whatDiseases_Treatment{i}; 

        [RANDOMdrugs(:,j), geneNames] = give_randomDrug_null(diseaseName, disorderDrugs, allDrugs, whatSelection);

    end
    RANDOMdrugs_treatment{i} = RANDOMdrugs;
end

switch whatSelection
    case 'random' % these are for old version of the data
        fileName = sprintf('DataOutput_2024/nulls_%d_%stargets_randomDrugR.mat', numNulls,params.whatDrugTargets); 
    case 'proportional' % these are for old version of the data
        fileName = sprintf('DataOutput_2024/nulls_%d_%stargets_randomDrugP.mat', numNulls,params.whatDrugTargets); 
    case 'proportionalPsych' % these are used for drugbank data
        fileName = sprintf('DataOutput_2024/nulls_%d_%stargets_randomDrugP_%s_drugbank_psych.mat', numNulls,params.whatDrugTargets, params.whatTargets); 
    case 'drugbank'  % these are used for drugbank data
        fileName = sprintf('DataOutput_2024/nulls_%d_%stargets_randomDrugR_%s_drugbank.mat', numNulls,params.whatDrugTargets, params.whatTargets); 
end

save(fileName, 'RANDOMdrugs_treatment', 'whatDiseases_Treatment', 'geneNames', 'params');
end