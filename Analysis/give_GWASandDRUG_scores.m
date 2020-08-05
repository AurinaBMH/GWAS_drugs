function [geneWeightsGWAS_ALL, drugScores_ord] = give_GWASandDRUG_scores(whatGWAS)

whatThreshold = 'BF'; 
similarityTypes = {'Adult_brain', 'AllenMeanCoexpMapped', 'AllenMeanCoexpeQTLbrain', ...
        'Astro', 'Fetal_brain', 'MAGMAdefault', 'Neuro', ...
        'PPI_eQTLbrain_th0', 'PPI_eQTLbrain_th400', 'PPI_eQTLbrain_th600', 'PPI_eQTLbrain_th900', ...
        'PPI_mapped_th0', 'PPI_mapped_th400', 'PPI_mapped_th600', 'PPI_mapped_th900', ...
        'eQTLHeart_Left_Ventricle', 'eQTLLiver', 'eQTLWhole_Blood', 'eQTLbrain'};
PPImeasures_names = {'numPPIneighbors1','percPPIneighbors1','weiPPIneighbors1','expWeiPPIneighbors1', 'numPPIneighbors2','percPPIneighbors2', 'weiPPIneighbors2','expWeiPPIneighbors2'};    

whatDiseases_GWAS = {'ADHD', 'MDD2', 'SCZ', 'BIP2', 'DIABETES', 'HF'};
whatDiseases_Treatment = {'ADHD','BIP','SCZ','MDD','cardiology','diabetes'};


Dname = whatDiseases_GWAS{1};
if strcmp(Dname, 'HF')
    Dname = 'cardiology';
end
Dname = Dname(isstrprop(Dname,'alpha'));
takeVal = contains(whatDiseases_Treatment, Dname, 'IgnoreCase',true);


%-------------------------------------------------------------------------------
% Load treatment weights of each gene implicated in each disorder:
[geneNamesDrug,drugScoresAll] = GiveMeNormalizedScoreVectors(whatDiseases_Treatment,'Drug');

numScores = length(find(contains(similarityTypes,'PPI')))*length(PPImeasures_names)+length(find(~contains(similarityTypes,'PPI')));
geneWeightsGWAS = cell(numScores,1);
k=1;
for s=1:length(similarityTypes)
    if contains(similarityTypes{s},'PPI')
        for p=1:length(PPImeasures_names)
            
            whatProperty = PPImeasures_names{p};
            [geneNamesGWAS,geneWeightsGWAS{k}] = GiveMeNormalizedScoreVector(whatGWAS,'GWAS',similarityTypes{s},whatProperty, whatThreshold);
            k=k+1;
        end
    else
        if ~contains(similarityTypes{s},'Allen')
            whatProperty = 'P';
        elseif contains(similarityTypes{s},'Allen')
            whatProperty = 'r';
        end
        [geneNamesGWAS,geneWeightsGWAS{k}] = GiveMeNormalizedScoreVector(whatGWAS,'GWAS',similarityTypes{s},whatProperty, whatThreshold);
        k=k+1;
    end
end

% Combine two datasets on overlap:
[geneNames,ia,ib] = intersect(geneNamesGWAS,geneNamesDrug,'stable');
geneWeightsGWAS_ord = cell(length(geneWeightsGWAS),1); 
for k=1:length(geneWeightsGWAS)
    
    geneWeightsGWAS_ord{k} = geneWeightsGWAS{k}(ia,:);
    
end
drugScores_ord = drugScoresAll(ib,:);
drugScores_ord = drugScores_ord(:,takeVal);

geneWeightsGWAS_ALL = horzcat(geneWeightsGWAS_ord{:}); 

% make a vector of names for GWAS measures; 

end

