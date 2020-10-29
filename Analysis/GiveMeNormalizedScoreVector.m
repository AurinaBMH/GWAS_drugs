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

params = SetDefaultParams();
whatNorm = params.whatNorm; 

% this norm here doesn't matter as much, when computing the dot-product,
% scores are re-normalized, so that's where it matters
% but lley's use consistent
%-------------------------------------------------------------------------------

switch whatMeasurement
case 'Drug'
    normalizeWithinDrugs = true;
    indicatorTable = ImportTreatmentLists(normalizeWithinDrugs, params.whatDrugTargets, params.whatTargets);
    geneNames = indicatorTable.Row;
    geneWeights = indicatorTable.(whatDisease);
case 'GWAS'
    % Load data:
    fileName = sprintf('resultsTable_%s_%s_%s_%s_drugbank.mat',whatDisease, whatThreshold, params.whatDrugTargets, params.whatTargets);
    load(fileName,'geneScores');
    fprintf(1,'Loaded gene scores for %s from %s\n',whatDisease,fileName);
    geneNames = geneScores.gene;
    
    % try and keep only thresholded set of genes, give 0 to all others
%     load('GWAS_disordersMAGMA.mat')
%     listGENESmapped = DISORDERlist.MAGMAdefault.(whatDisease);
%     pThr_m = 0.05/size(listGENESmapped,1); % Bonf correction for the number of genes in the list
%     allMappedDiseaseGenes = listGENESmapped.GENENAME(listGENESmapped.P<pThr_m);
%     
%     isGWAS = ismember(geneNames,allMappedDiseaseGenes);

    if strcmp(similarityType, 'MAGMAdefault') || strcmp(similarityType, 'eQTLbrain')
        % BF correction over 2000 genes leaves almost no genes, use p<0.01
        % threshold, this is -log10(p)>2; 
        if strcmp(whatProperty, 'P')
            geneWeights = geneScores.(similarityType).(whatProperty);
            %isGWAS = geneScores.(similarityType).(whatProperty)>-log10(0.05);
            % null all scores with p>0.05 and keep weights for others; 
            notGWAS = geneScores.(similarityType).(whatProperty)<-log10(0.05);
            %isGWAS = (10.^-(geneScores.(similarityType).(whatProperty)))<0.05/length(geneNames); 
            geneWeights(notGWAS) = 0;
            %geneWeights = double(isGWAS); 
        else
            geneWeights = geneScores.(similarityType).(whatProperty);
        end
    elseif contains(similarityType, 'Allen')
        if strcmp(whatProperty, 'zval')
            % take absolute of z-value
            % geneWeights = -geneScores.(similarityType).(whatProperty);
            geneWeights = abs(geneScores.(similarityType).(whatProperty));
        else
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
