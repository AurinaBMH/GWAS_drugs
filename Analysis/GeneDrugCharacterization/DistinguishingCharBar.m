function DistinguishingCharBar(similarityType,whatProperty)

if nargin < 1
    similarityType = 'DNA';
end
if nargin < 2
    whatProperty = 'propGWAS';
    % 'numGWASMapped','numLDSNPs','percPPIneigh1Mapped','percPPIneigh1LD','AllenMeanCoexpMapped','AllenMeanCoexpLD'
end

whatDiseases_GWAS = {'ADHD','BIP','SZP','MDD','diabetes'};
whatDiseases_Treatment = {'ADHD','BIP','SZP','MDD','pulmonary','cardiology','gastro','diabetes'};

%-------------------------------------------------------------------------------
% Load in default parameters:
params = SetDefaultParams();
whatScore = params.whatScore;

%-------------------------------------------------------------------------------
numDiseases_Treatment = length(whatDiseases_Treatment);
numDiseases_GWAS = length(whatDiseases_GWAS);

%-------------------------------------------------------------------------------
% Load treatment weights of each gene implicated in each disorder:
normalizeWithinDrugs = true; % weight genes lower if they occur in drugs with large numbers of gene targets
[indicatorTable,percIndicatorTable] = ImportTreatmentLists(normalizeWithinDrugs);

%===============================================================================
f = figure('color','w');
for i = 1:numDiseases_GWAS
    % Load data:
    fileName = sprintf('resultsTable_%s.mat',whatDiseases_GWAS{i});
    load(fileName,'geneScores');
    fprintf(1,'Loaded gene scores from %s\n',fileName);
    switch similarityType
    case 'DNA'
        % numSNPs, numGWAS, percGWAS, numLD, percLD
        geneWeights_GWAS = geneScores.DNA.(whatProperty);
    case 'PPI_weighted'
        % numPPIneighbors1
        % percPPIneighbors1
        % weiPPIneighbors1
        % medianPPIDistance
        % meanPPIDistance
        geneWeights_GWAS = geneScores.PPI_mapped_weighted.(whatProperty);
    case 'PPI_th0'
        geneWeights_GWAS = geneScores.PPI_mapped_th0.(whatProperty);
    case 'PPI_th400'
        geneWeights_GWAS = geneScores.PPI_mapped_th400.(whatProperty);
    case 'PPI_th600'
        geneWeights_GWAS = geneScores.PPI_mapped_th600.(whatProperty);
    case 'PPI_th900'
        geneWeights_GWAS = geneScores.PPI_mapped_th900.(whatProperty);
    case 'Expression'
        geneWeights_GWAS = geneScores.AllenMeanCoexpMapped;
    otherwise
        error('Unknown similarity type: ''%s''',similarityType);
    end

    % Combine two datasets on overlap:
    [~,ia,ib] = intersect(geneScores.gene,indicatorTable.Properties.RowNames,'stable');
    geneWeights_GWAS = geneWeights_GWAS(ia);
    indicatorTable = indicatorTable(ib,:);
    % resultsTable = resultsTable(ia,:);

    fprintf(1,'%u matching (/%u %s GWAS); (/%u with treatment weights)\n',...
        length(ia),length(geneWeights_GWAS),whatDiseases_GWAS{i},height(indicatorTable));

    %-------------------------------------------------------------------------------
    % Get scores for the property of interest:
    rhos = zeros(numDiseases_Treatment,1);
    isSig = zeros(numDiseases_Treatment,1);
    switch whatScore
    case 'Kendall'
        addNull = false;
        for k = 1:numDiseases_Treatment
            geneWeights_treatment = indicatorTable.(whatDiseases_Treatment{k});
            [rhos(k),p] = corr(geneWeights_GWAS,geneWeights_treatment,...
                                    'type','Kendall','rows','pairwise');
            isSig(k) = (p < 0.05);
        end
    case 'weightedSum'
        addNull = true;
        for k = 1:numDiseases_Treatment
            theDisease = whatDiseases_Treatment{k};
            geneWeights_treatment = indicatorTable.(theDisease);
            rhos(k) = ComputeWeightedSum(geneWeights_treatment,geneWeights_GWAS,theDisease,false);
        end

        % Generate null distributions:
        numNulls = 1000;
        nullScores = zeros(numNulls,1);
        whatNull = 'randomWeight'; % randomWeight, randomDisease
        for k = 1:numNulls
            switch whatNull
            case 'randomWeight'
                geneWeightsRand = rand(height(indicatorTable),1);
                nullScores(k) = ComputeWeightedSum(geneWeightsRand,geneWeights_GWAS,theDisease,false);
            case 'randomDisease'
                % Shuffle weights taken from a random disease (pooled nulls):
                % [could be done individually for each particular weighting if needed]
                theDiseaseInd = randi(numDiseases_Treatment,1);
                theDisease = whatDiseases_Treatment{theDiseaseInd};
                geneWeightsRand = indicatorTable.(theDisease);
                nullScores(k) = ComputeWeightedSum(geneWeightsRand,geneWeights_GWAS,theDisease,true);
            end
        end
        % Compute p-values:
        for k = 1:numDiseases_Treatment
            isSig(k) = (mean(rhos(k) < nullScores) < 0.05);
        end
    end
    [rhos,ix] = sort(rhos,'descend');

    %---------------------------------------------------------------------------
    ax = subplot(1,5,i); hold on
    b = bar(rhos);
    if addNull && ~all(isnan(nullScores))
        % Add null distribution:
        [ff,x] = ksdensity(nullScores,linspace(min(nullScores),max(nullScores),500),'function','pdf');
        ff = 0.8*ff/max(ff);
        plot(ones(2,1)*(numDiseases_Treatment+1)+ff,x,'k');
        plot(ones(2,1)*(numDiseases_Treatment+1)-ff,x,'k');
        ax.XTick = 1:numDiseases_Treatment+1;
        ax.XTickLabel = {whatDiseases_Treatment{ix},'null'};
        % Add horizontal lines to aid comparison to null:
        null_p50 = quantile(nullScores,0.5);
        plot([1,numDiseases_Treatment+1],ones(2,1)*null_p50,':k')
        null_p95 = quantile(nullScores,0.95);
        plot([1,numDiseases_Treatment+1],ones(2,1)*null_p95,'--k')
        if range(rhos) > 0
            maxLim = max(rhos)*1.1;
            if maxLim < null_p95
                maxLim = null_p95*1.1;
            end
            minLim = min(rhos)*0.9;
            null_p10 = quantile(nullScores,0.1);
            if minLim > null_p10;
                minLim = null_p10*0.9;
            end
            ax.YLim = [minLim,maxLim];
        end
    else
        ax.XTick = 1:numDiseases_Treatment;
        ax.XTickLabel = whatDiseases_Treatment(ix);
    end
    ax.XTickLabelRotation = 45;
    xlabel('Disease treatment')
    ylabel(sprintf('%s similarity',whatScore))
    title(sprintf('%s-%s',whatDiseases_GWAS{i},whatProperty),'interpreter','none')
    cMapGeneric = BF_getcmap('set2',numDiseases_Treatment,false);
    for k = find(~isSig)'
        cMapGeneric(k,:) = brighten(cMapGeneric(k,:),0.8);
    end
    b.CData = cMapGeneric(ix,:);
    b.FaceColor = 'flat';
end

f.Position = [471,945,1830,318];

%-------------------------------------------------------------------------------
function rho = ComputeWeightedSum(v1,v2,theName,doShuffle)
    % Compute the weighted sum of two vectors, v1 & v2
    if nargin < 3
        theName = '';
    end
    if nargin < 4
        doShuffle = false;
    end
    minGoodProp = 0.4; % Require fewer than 40% bad values after matching indices

    % Check for and filter bad values:
    r = ~isnan(v1) & ~isnan(v2);
    if mean(r) > minGoodProp
        warning('More than 40% bad values when comparing %s',theName)
    end
    v1_good = v1(r);
    v2_good = v2(r);
    if doShuffle
        % Shuffle weights on v1:
        v1_good = v1_good(randperm(length(v1_good)));
    end

    % Normalize v1 and v2 to sum to 1:
    if any(v1_good < 0) || any(v2_good < 0)
        warning('Weight vectors contain negative values')
    end
    if all(v1_good > 0)
        v1_good = v1_good/norm(v1_good,2);
    end
    if all(v2_good > 0)
        v2_good = v2_good/norm(v2_good,2);
    end

    % Compute the score:
    rho = sum(v1_good.*v2_good);
end

end
