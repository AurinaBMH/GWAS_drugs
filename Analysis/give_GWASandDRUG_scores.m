function [geneWeightsGWAS_ALL, drugScores_ord, similarityTypes, PPImeasures_names, measureNames] = give_GWASandDRUG_scores(whatGWAS, whatMeasures)

whatThreshold = 'BF'; 
params = SetDefaultParams();

switch whatMeasures
    case 'all'
similarityTypes = [{'AllenMeanCoexpMapped'}, {'AllenMeanCoexpeQTLbrain'}, ...
        {'PPI_eQTLbrain_th0'}, {'PPI_eQTLbrain_th400'}, {'PPI_eQTLbrain_th600'}, {'PPI_eQTLbrain_th900'}, ...
        {'PPI_mapped_th0'}, {'PPI_mapped_th400'}, {'PPI_mapped_th600'}, {'PPI_mapped_th900'}, ...
        params.whatANNOT(:)'];
PPImeasures_names = {'numPPIneighbors1','percPPIneighbors1'};  

    case 'reduced'
similarityTypes = {'Adult_brain', 'MAGMAdefault', 'eQTLbrain','AllenMeanCoexpMapped', ...
        'PPI_mapped_th900', 'PPI_eQTLbrain_th400','eQTLLiver', 'eQTLPancreas', 'eQTLColon_Transverse', 'eQTLHeart_Left_Ventricle'}';
PPImeasures_names = {'percPPIneighbors1'};   
end
    
 
whatDiseases_GWAS =params.whatGWAS; 
whatDiseases_Treatment = params.whatDiseases_Treatment; 

Dname = whatDiseases_GWAS{1};
Dname = Dname(isstrprop(Dname,'alpha'));
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
            whatProperty = 'r';
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

% make a vector of names for GWAS measures; 

end

