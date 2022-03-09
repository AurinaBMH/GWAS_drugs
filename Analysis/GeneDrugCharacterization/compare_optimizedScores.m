function [Ptable, measureTable] = compare_optimizedScores(whatDiseases_GWAS, whatMeasures, whatNull)
%X - predictor data: all gene scores from GWAS
%y - response: drug score
if nargin < 1
    params = SetDefaultParams();
    whatDiseases_GWAS = params.whatGWAS_2022;
end
if nargin < 2
    whatMeasures = 'allPsych'; % 'reduced' or 'all';
end

if nargin < 3
    whatNull = 'randomDrugR_all_drugbank';
end

params = SetDefaultParams();
whatNorm = params.whatNorm;

numGWAS = length(whatDiseases_GWAS);

for i = 1:numGWAS
    
    whatGWAS = whatDiseases_GWAS{i};
    [geneWeightsGWAS_ALL, drugScores_ord, similarityTypes, PPImeasures_names] = give_GWASandDRUG_scores(whatGWAS, whatMeasures);
    
    % remove rows that are all NaN and measure names
    INDnan = all(isnan(geneWeightsGWAS_ALL),1);
    geneWeightsGWAS_ALL(:, INDnan) = [];
    
    % plot randomNullP-based p-values for each mapping method
    Dname = whatGWAS(isstrprop(whatGWAS,'alpha'));
    
    [~, diseaseResultsP, ~,~, measureNames] = compareGWASvsDRUGmatches({whatGWAS}, whatNull, Dname, similarityTypes, PPImeasures_names, whatMeasures);
    close all;
    
    diseaseResultsP(diseaseResultsP==0) = 1/params.numNull;
    Pvals_mapp = -log10(diseaseResultsP(~isnan(diseaseResultsP(:))))';

    % apply linear regression
    mdl = fitlm(geneWeightsGWAS_ALL,drugScores_ord);
    ypred = predict(mdl,geneWeightsGWAS_ALL);
    ypredNorm = normalizeScoreVector(ypred, whatNorm);
    
    Pval_comb = compare_to_null(whatGWAS, ypredNorm, drugScores_ord, geneWeightsGWAS_ALL, whatNull, whatNorm);
    
    if Pval_comb==0
        Pval_comb = 1/params.numNull;
    end
    pPLOT = -log10(Pval_comb);
    
    measureTable = [measureNames; 'Combined'];
    PvalsTable = [Pvals_mapp, pPLOT];
    
    Ptable.(whatGWAS).measureName = measureTable; 
    Ptable.(whatGWAS).Pvals = PvalsTable; 

    
end

end


