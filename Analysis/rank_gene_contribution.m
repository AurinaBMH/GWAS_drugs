% this function ranks each gene by its contribution to the matching between
% GWAS and drug target
function [Tp, Tdot] = rank_gene_contribution(whatGWAS, whatDrugs, similarityType, whatNull)

if nargin < 4
    whatNull = 'randomDrugR_all_drugbank'; 
end

% load null data
fileName = sprintf('nulls_5000_2020targets_%s.mat', whatNull); 
load(fileName, 'whatDiseases_Treatment', 'RANDOMdrugs_treatment', 'geneNames');


params = SetDefaultParams();
whatDiseases_Treatment_ALL = params.whatDiseases_Treatment_ALL;

% get drug scores and select treatments
[geneNamesDrug,drugScoresAll] = GiveMeNormalizedScoreVectors(whatDiseases_Treatment_ALL,'Drug');
[~, INDtreatment] = intersect(whatDiseases_Treatment_ALL, whatDrugs);

% get GWAS scores
if contains(similarityType,'PPI')
    whatProperty = 'percPPIneighbors1';
else
    if ~contains(similarityType,'Allen')
        whatProperty = 'P';
    elseif contains(similarityType,'Allen')
        whatProperty = 'r';
    end
end
[geneNamesGWAS,geneWeightsGWAS] = GiveMeNormalizedScoreVector(whatGWAS,'GWAS',similarityType, whatProperty, params.whatThreshold);

% Combine two datasets on overlap:
[geneNames_matched,ia,ib] = intersect(geneNamesGWAS,geneNamesDrug);
geneWeightsGWAS = geneWeightsGWAS(ia);
geneWeightsDRUG = drugScoresAll(ib,INDtreatment);

% Combine nulls and real data on overlap
[~,~,ir] = intersect(geneNames_matched, geneNames, 'stable');
[~, INDnull] = intersect(whatDiseases_Treatment, whatDrugs);
geneWeightsNull = RANDOMdrugs_treatment{INDnull}(ir,:); 

numGenes = length(geneWeightsGWAS); 
numNulls = size(geneWeightsNull,2); 

% for each gene, get the dot product in real data
% get dot product for null
% compare real vs null, get p-value
% save p-value, rank based on p-values

% for a different GWAS on the same drugs, everything is rescaled by GWAS
% values, but the ranking is the same
dotP = nan(numGenes,1); 
dotP_null = nan(numGenes, numNulls); 
p_val = nan(numGenes,1); 
for g=1:numGenes
    
    dotP(g) = geneWeightsGWAS(g).*geneWeightsDRUG(g);
    % if the mathch in real data is 0, it's not contributing
    if dotP(g)~=0 && ~isnan(dotP(g))
        % get pvalue for others
        dotP_null(g,:) = geneWeightsGWAS(g).*geneWeightsNull(g,:);
        p_val(g) = mean(dotP_null(g,:)>dotP(g));
    end
end
[~, iP] = sort(p_val, 'ascend'); 
% rank p_vals and give out a ranked table of all genes
Tp = table; 
Tp.Gene = geneNames(iP); 
Tp.Pval = p_val(iP); 

[~, iD] = sort(dotP, 'descend', 'MissingPlacement','last'); 
Tdot = table; 
Tdot.Gene = geneNames(iD); 
Tdot.DotP = dotP(iD); 

end
