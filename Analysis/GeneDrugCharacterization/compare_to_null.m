function pVals = compare_to_null(whatGWAS, geneWeightGWAS, drugScores_DIS, geneWeightsGWAS_all, whatNull, whatNorm)

if nargin < 5
    whatNull = 'randomDrugR_all_drugbank';
    whatNorm = 2; 
end

params = SetDefaultParams();
whatGWAS = whatGWAS(isstrprop(whatGWAS,'alpha'));
load(sprintf('nulls_5000_%stargets_%s.mat', params.whatDrugTargets, whatNull), 'RANDOMdrugs_treatment', 'whatDiseases_Treatment', 'geneNames');


% choose null for the corresponding disorder
[~, l] = intersect(whatDiseases_Treatment,whatGWAS); 
% Generate null distributions:
numNulls = params.numNull; 
rhos = ComputeDotProduct(drugScores_DIS,geneWeightGWAS);
nullScores = zeros(numNulls,1);

for k = 1:numNulls
    % separate set of nulls for each drug target list
    if contains(whatNull, 'randomTarget')
        % is the actual match higher than a match with random gene score assignment
        % Shuffle weights taken from each drug list individually
        nullScores(k) = ComputeDotProduct(drugScores_DIS,geneWeightGWAS, true);
        % randomise v1 within ComputeDotProduct
%     elseif contains(whatNull, 'randomDrug')
%         % for each disease get a random set of drugs that is the same size as
%         drugScores_DISrand = RANDOMdrugs_treatment{l}(:,k);
%         nullScores(k) = ComputeDotProduct(drugScores_DISrand,geneWeightGWAS);
    elseif contains(whatNull, 'randomDrug')
        % select null drug vector
        drugScores_DISrand = RANDOMdrugs_treatment{l}(:,k);
        
        % apply linear regression based on all measures
        mdl = fitlm(geneWeightsGWAS_all,drugScores_DISrand);
        ypred = predict(mdl,geneWeightsGWAS_all);
        ypredNorm = normalizeScoreVector(ypred, whatNorm);
        
        % get the score
        nullScores(k) = ComputeDotProduct(drugScores_DISrand,ypredNorm);

    end
end

if ~isnan(rhos)
    % Compute p-values: based on separate nulls
    pVals = mean(rhos < nullScores);
else
    pVals = NaN;
end
    

end


