function [geneWeightsGWAS_ALL, drugScores_ord, similarityTypes, PPImeasures_names, measureNames, whatDiseases_Treatment] = give_GWASandDRUG_scores(whatGWAS, whatMeasures)

whatThreshold = 'BF'; 
params = SetDefaultParams();

switch whatMeasures
    case 'allPsych'
        similarityTypes = params.whatANNOT_psych;
        PPImeasures_names = {'numPPIneighbors1','percPPIneighbors1'};
        whatDiseases_Treatment = params.whatDiseases_Treatment; 
        
    case 'allBody'
        similarityTypes = params.whatANNOT_body;
        PPImeasures_names = {'numPPIneighbors1','percPPIneighbors1'};
        whatDiseases_Treatment = params.whatDiseases_Treatment_body; 
        
    case 'reduced'
        similarityTypes = params.whatANNOT_reduced; 
        PPImeasures_names = {'percPPIneighbors1'};
        whatDiseases_Treatment = params.whatDiseases_Treatment; 
        
    case 'all'
        similarityTypes = params.whatANNOT_all; 
        PPImeasures_names = {'numPPIneighbors1','percPPIneighbors1'};
        whatDiseases_Treatment = params.whatDiseases_Treatment_label_ALL; 
        
    case 'reduced5'
        similarityTypes = params.whatANNOT_reduced5; 
        PPImeasures_names = {'percPPIneighbors1'};
        whatDiseases_Treatment = params.whatDiseases_Treatment; 
        
end

Dname = whatGWAS(isstrprop(whatGWAS,'alpha'));
takeVal = contains(whatDiseases_Treatment, Dname, 'IgnoreCase',true);
%-------------------------------------------------------------------------------
% Load treatment weights of each gene implicated in each disorder:
[geneNamesDrug,drugScoresAll] = GiveMeNormalizedScoreVectors(whatDiseases_Treatment,'Drug');

numScores = length(find(contains(similarityTypes,'PPI')))*length(PPImeasures_names)+length(find(~contains(similarityTypes,'PPI')));
geneWeightsGWAS = cell(numScores,1);
measureNames = cell(numScores,1);
k=1;
for s=1:length(similarityTypes)
    if contains(similarityTypes{s},'PPI')
        for p=1:length(PPImeasures_names)
            
            whatProperty = PPImeasures_names{p};
            [geneNamesGWAS,geneWeightsGWAS{k}] = GiveMeNormalizedScoreVector(whatGWAS,'GWAS',similarityTypes{s},whatProperty, whatThreshold);
            measureNames{k} = [similarityTypes{s}, whatProperty]; 
            k=k+1;
        end
    else
        if ~contains(similarityTypes{s},'Allen')
            whatProperty = 'P';
        elseif contains(similarityTypes{s},'Allen')
            whatProperty = 'zval';
        end
        [geneNamesGWAS,geneWeightsGWAS{k}] = GiveMeNormalizedScoreVector(whatGWAS,'GWAS',similarityTypes{s},whatProperty, whatThreshold);
        measureNames{k} = [similarityTypes{s}, whatProperty]; 
        k=k+1;
    end
end

% Combine two datasets on overlap:
[~,ia,ib] = intersect(geneNamesGWAS,geneNamesDrug);

geneWeightsGWAS_ord = cell(length(geneWeightsGWAS),1); 
for k=1:length(geneWeightsGWAS)
    
    geneWeightsGWAS_ord{k} = geneWeightsGWAS{k}(ia,:);
    
end
drugScores_ord = drugScoresAll(ib,takeVal);

geneWeightsGWAS_ALL = horzcat(geneWeightsGWAS_ord{:}); 

end

