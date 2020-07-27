function [geneNames,geneWeightsNorm,geneWeights] = GiveMeNormalizedScoreVector(whatDisease,whatMeasurement,similarityType,whatProperty, whatThreshold)
%-------------------------------------------------------------------------------
% Inputs:
%-------------------------------------------------------------------------------
if nargin < 1
    whatDisease = 'SCZ';
end
if nargin < 2
    whatMeasurement = 'GWAS';
end
if nargin < 3
    similarityType = 'MAGMAdefault';
    fprintf('using %s similarity type by DEFAULT\n', similarityType)
end
if nargin < 4
    if contains(similarityType, 'PPI')
        whatProperty = 'percPPIneighbors1';
    else
        whatProperty = 'P';
    end
end
if nargin <5
     whatThreshold = 'BF'; 
     fprintf('using %s threshold by DEFAULT\n', whatThreshold)
end
whatNorm = 2;
% this norm here doesn't matter as much, when computing the dot-product,
% scores are re-normalized, so that's where it matters
% but lley's use consistent
%-------------------------------------------------------------------------------

switch whatMeasurement
case 'Drug'
    normalizeWithinDrugs = true;
    indicatorTable = ImportTreatmentLists(normalizeWithinDrugs);
    geneNames = indicatorTable.Row;
    geneWeights = indicatorTable.(whatDisease);
case 'GWAS'
    % Load data:
    fileName = sprintf('resultsTable_%s_%s.mat',whatDisease, whatThreshold);
    load(fileName,'geneScores');
    fprintf(1,'Loaded gene scores for %s from %s\n',whatDisease,fileName);
    geneNames = geneScores.gene;

    
    if strcmp(similarityType, 'AllenMeanCoexpMapped') || strcmp(similarityType, 'AllenMeanCoexpeQTLbrain')
        geneWeights = geneScores.(similarityType);
    else
        geneWeights = geneScores.(similarityType).(whatProperty);
    end
    
end

%-------------------------------------------------------------------------------
% Normalize (non-NaN elements) to unit vector as 2-norm:
geneWeightsNorm = normalizeScoreVector(geneWeights, whatNorm); 

end
