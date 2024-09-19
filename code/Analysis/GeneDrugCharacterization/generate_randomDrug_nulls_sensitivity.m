function generate_randomDrug_nulls_sensitivity(whatSelection)
if nargin<1
    whatSelection = 'drugbank'; 
end
% this function generates a set of random nulls by selecting the selected
% numnber of genes each time; Overall, the nulls for different numbers of
% drugs are very similar, could use only one null with the average number
% of drugs for comparison to all, but for now making a separate null for
% each drug list; 

numNulls = 10000; 
params = SetDefaultParams();
treatment_classes = params.whatDiseases_Treatment_classes;

[geneNamesGWAS,~] = GiveMeNormalizedScoreVector('ADHD3','GWAS','MAGMAdefault',params.geneScore, params.whatThreshold);
[~,~, disorderDrugs, allDrugs] = ImportTreatmentLists_sensitivity(true, 'sensitivity', params.whatTargets);

numGenes = length(geneNamesGWAS); 
RANDOMdrugs = zeros(numGenes, numNulls); 
RANDOMdrugs_treatment = cell(length(treatment_classes),1); 

for i = 1:length(treatment_classes)
    
    fprintf('Running nulls for %s\n', treatment_classes{i})
    
    for j = 1:numNulls
        diseaseName = treatment_classes{i}; 

        [RANDOMdrugs(:,j), geneNames] = give_randomDrug_null_sensitivity(diseaseName, disorderDrugs, allDrugs, whatSelection);

    end
    RANDOMdrugs_treatment{i} = RANDOMdrugs;
end

switch whatSelection
    case 'random' % these are for old version of the data
        fileName = sprintf('DataOutput_2024/nulls_%d_%stargets_randomDrugR_sensitivity.mat', numNulls,params.whatDrugTargets); 
    case 'proportional' % these are for old version of the data
        fileName = sprintf('DataOutput_2024/nulls_%d_%stargets_randomDrugP_sensitivity.mat', numNulls,params.whatDrugTargets); 
    case 'proportionalPsych' % these are used for drugbank data
        fileName = sprintf('DataOutput_2024/nulls_%d_%stargets_randomDrugP_%s_drugbank_psych_treatment_class.mat', numNulls,params.whatDrugTargets, params.whatTargets); 
    case 'drugbank'  % these are used for drugbank data
        fileName = sprintf('DataOutput_2024/nulls_%d_%stargets_randomDrugR_%s_drugbank_treatment_class.mat', numNulls,params.whatDrugTargets, params.whatTargets); 
end

save(fileName, 'RANDOMdrugs_treatment', 'treatment_classes', 'geneNames', 'params');
end