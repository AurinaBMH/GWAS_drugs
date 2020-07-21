function [geneNames,geneWeightsNorm,geneWeights] = GiveMeNormalizedScoreVector(whatDisease,whatMeasurement,similarityType,whatProperty, whatThreshold)
%-------------------------------------------------------------------------------
% Inputs:
%-------------------------------------------------------------------------------
if nargin < 1
    whatDisease = 'SZP';
end
if nargin < 2
    whatMeasurement = 'GWAS';
end
if nargin < 3
    similarityType = 'MAGMAdefault';
end
if nargin < 4
    if contains(similarityType, 'PPI')
        whatProperty = 'percPPIneighbors1';
    else
        whatProperty = 'ZSTAT';
    end
end
whatNorm = 1;

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
r = ~isnan(geneWeights);
geneWeightsNorm = geneWeights;
geneWeightsNorm(r) = geneWeightsNorm(r)/norm(geneWeights(r),whatNorm);

end
