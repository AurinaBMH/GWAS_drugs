function [geneNames,geneWeightsNorm,geneWeights] = GiveMeNormalizedScoreVector_sensitivity(whatDisease,whatMeasurement,similarityType,whatProperty, whatThreshold)
%---------------------------------------r----------------------------------------
% Inputs:
%-------------------------------------------------------------------------------
if nargin < 1
    whatDisease = 'BIP';
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

params = SetDefaultParams();
whatNorm = params.whatNorm; 

% this norm here doesn't matter as much, when computing the dot-product,
% scores are re-normalized, so that's where it matters
% but lley's use consistent
%-------------------------------------------------------------------------------

switch whatMeasurement
case 'Drug'
    normalizeWithinDrugs = true;
    indicatorTable = ImportTreatmentLists(normalizeWithinDrugs, 'sensitivity', params.whatTargets);
    geneNames = indicatorTable.Row;
    geneWeights = indicatorTable.(whatDisease);
case 'GWAS'
    % Load data:
    fileName = sprintf('DataOutput_2024/resultsTable_%s_%s_%s_%s_drugbank.mat',whatDisease, whatThreshold, params.whatDrugTargets, params.whatTargets);
    load(fileName,'geneScores');
    fprintf(1,'Loaded gene scores for %s from %s\n',whatDisease,fileName);
    geneNames = geneScores.gene;
    
    if strcmp(similarityType, 'MAGMAdefault') || strcmp(similarityType, 'eQTLbrain')
        if strcmp(whatProperty, 'P') && params.DOthreshold
            % if there is a selection to threshold, remove scores for genes
            % that had p>0.05, keep for others;
            geneWeights = geneScores.(similarityType).(whatProperty);
            notGWAS = geneScores.(similarityType).(whatProperty)<-log10(0.05);
            geneWeights(notGWAS) = 0;
        else
            geneWeights = geneScores.(similarityType).(whatProperty);
        end
    elseif contains(similarityType, 'Allen')
        if strcmp(whatProperty, 'zval') && strcmp(params.AHBAmeasure, '-z')
            % take negative of z
            geneWeights = -geneScores.(similarityType).(whatProperty);
        elseif strcmp(whatProperty, 'zval') && strcmp(params.AHBAmeasure, 'abs(z)')
            % take abs of z
            geneWeights = abs(geneScores.(similarityType).(whatProperty));
        else
            % take z
            geneWeights = geneScores.(similarityType).(whatProperty);
        end
    else
        
        geneWeights = geneScores.(similarityType).(whatProperty);
    end

end

%-------------------------------------------------------------------------------
% Normalize (non-NaN elements) to unit vector as 2-norm:
geneWeightsNorm = normalizeScoreVector(geneWeights, whatNorm);

end
