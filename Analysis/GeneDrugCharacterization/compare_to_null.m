function pVals = compare_to_null(whatGWAS, geneWeightsGWAS, drugScores_DIS, whatNull)

if nargin < 4
    whatNull = 'randomDrugP';
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
        switch whatNull
            
            case 'randomTarget' % is the actual match higher than a match with random gene score assignment
                % Shuffle weights taken from each drug list individually
                nullScores(k) = ComputeDotProduct(drugScores_DIS,geneWeightsGWAS, true);
                % randomise v1 within ComputeDotProduct
            case {'randomDrugP','randomDrugR'}  % for each disease get a random set of drugs that is the same size as
                drugScores_DISrand = RANDOMdrugs_treatment{l}(:,k);
                nullScores(k) = ComputeDotProduct(drugScores_DISrand,geneWeightsGWAS);
                
        end
    end
    % Compute p-values: based on separate nulls
    isSig = (mean(rhos < nullScores) < 0.05);
    pVals = mean(rhos < nullScores);

end


