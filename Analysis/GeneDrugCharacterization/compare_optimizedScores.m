function [coef_lasso, coef_linear, sigMeasures] = compare_optimizedScores(whatMeasures)
%X - predictor data: all gene scores from GWAS
%y - response: drug score
if nargin < 1
    whatMeasures = 'reduced'; % 'reduced' or 'all'; 
end

params = SetDefaultParams();
whatDiseases_GWAS = params.whatGWAS; 
whatNorm = params.whatNorm; 

numGWAS = length(whatDiseases_GWAS);

% run lasso for each GWAS
for i = 1:numGWAS
    
    whatGWAS = whatDiseases_GWAS{i};
    [geneWeightsGWAS_ALL, drugScores_ord, measureNames] = give_GWASandDRUG_scores(whatGWAS, whatMeasures);
    
    % remove rows that are all NaN and measure names
    INDnan = find(all(isnan(geneWeightsGWAS_ALL),1));
    geneWeightsGWAS_ALL(:, INDnan) = [];
    measureNames(INDnan) = [];
    
    % plot randomNullP-based p-values for each mapping method
    
    
    % apply linear regression
    mdl = fitlm(geneWeightsGWAS_ALL,drugScores_ord);
    ypred = predict(mdl,geneWeightsGWAS_ALL); 
    
    ypredNorm = normalizeScoreVector(ypred, whatNorm);
    
    % plot randomNullP-based p-value for combined measure
    
    
    
end

end
    
    
    