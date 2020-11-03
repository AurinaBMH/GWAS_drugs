function pVals = compare_to_null(whatGWAS, geneWeightsGWAS, drugScores_DIS, whatNull)

if nargin < 4
    whatNull = 'randomDrugR_all_drugbank';
end

params = SetDefaultParams();
whatGWAS = whatGWAS(isstrprop(whatGWAS,'alpha'));
load(sprintf('nulls_5000_%stargets_%s.mat', params.whatDrugTargets, whatNull), 'RANDOMdrugs_treatment', 'whatDiseases_Treatment', 'geneNames');


% choose null for the corresponding disorder
[~, l] = intersect(whatDiseases_Treatment,whatGWAS); 
% Generate null distributions:
numNulls = params.numNull; 
rhos = ComputeDotProduct(drugScores_DIS,geneWeightsGWAS);
nullScores = zeros(numNulls,1);

for k = 1:numNulls
    % separate set of nulls for each drug target list
    if contains(whatNull, 'randomTarget')
        % is the actual match higher than a match with random gene score assignment
        % Shuffle weights taken from each drug list individually
        nullScores(k) = ComputeDotProduct(drugScores_DIS,geneWeightsGWAS, true);
        % randomise v1 within ComputeDotProduct
    elseif contains(whatNull, 'randomDrug')
        % for each disease get a random set of drugs that is the same size as
        drugScores_DISrand = RANDOMdrugs_treatment{l}(:,k);
        nullScores(k) = ComputeDotProduct(drugScores_DISrand,geneWeightsGWAS);
        
    end
end
if ~isnan(rhos)
    % Compute p-values: based on separate nulls
    pVals = mean(rhos < nullScores);
else
    pVals = NaN;
end
    

end


