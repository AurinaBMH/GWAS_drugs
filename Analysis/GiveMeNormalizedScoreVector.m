function [geneNames,geneWeightsNorm,geneWeights] = GiveMeNormalizedScoreVector(whatDisease,whatMeasurement,similarityType,whatProperty)
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
    similarityType = 'DNA';
end
if nargin < 4
    whatProperty = 'percGWAS';
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
    fileName = sprintf('resultsTable_%s.mat',whatDisease);
    load(fileName,'geneScores');
    fprintf(1,'Loaded gene scores for %s from %s\n',whatDisease,fileName);
    geneNames = geneScores.gene;
    switch similarityType
    case 'DNA'
        % numSNPs, numGWAS, percGWAS, numLD, percLD
        geneWeights = geneScores.DNA.(whatProperty);
    case 'PPI_weighted'
        % numPPIneighbors1
        % percPPIneighbors1
        % weiPPIneighbors1
        % medianPPIDistance
        % meanPPIDistance
        geneWeights = geneScores.PPI_mapped_weighted.(whatProperty);
    case 'PPI_th0'
        geneWeights = geneScores.PPI_mapped_th0.(whatProperty);
    case 'PPI_th400'
        geneWeights = geneScores.PPI_mapped_th400.(whatProperty);
    case 'PPI_th600'
        geneWeights = geneScores.PPI_mapped_th600.(whatProperty);
    case 'PPI_th900'
        geneWeights = geneScores.PPI_mapped_th900.(whatProperty);
    case 'Expression'
        geneWeights = geneScores.AllenMeanCoexpMapped;
    otherwise
        error('Unknown similarity type: ''%s''',similarityType);
    end
end

%-------------------------------------------------------------------------------
% Normalize (non-NaN elements) to unit vector as 2-norm:
r = ~isnan(geneWeights);
geneWeightsNorm = geneWeights;
geneWeightsNorm(r) = geneWeightsNorm(r)/norm(geneWeights(r),whatNorm);

end
