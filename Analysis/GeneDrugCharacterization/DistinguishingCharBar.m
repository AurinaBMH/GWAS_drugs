function DistinguishingCharBar(similarityType,whatProperty)

if nargin < 1
    similarityType = 'DNA';
end
if nargin < 2
    whatProperty = 'percGWAS';
    % 'numGWASMapped','numLDSNPs','percPPIneigh1Mapped','percPPIneigh1LD','AllenMeanCoexpMapped','AllenMeanCoexpLD'
end
addNull = true;

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
[geneNamesDrug,drugScoresAll] = GiveMeNormalizedScoreVectors(whatDiseases_Treatment,'Drug');
numDrugScores = length(drugScoresAll);

%===============================================================================
f = figure('color','w');
ax = cell(numDiseases_GWAS,1);
for i = 1:numDiseases_GWAS
    whatDisease = whatDiseases_GWAS{i};
    [geneNamesGWAS,geneWeightsGWAS] = GiveMeNormalizedScoreVector(whatDisease,'GWAS',similarityType,whatProperty);

    % Combine two datasets on overlap:
    [geneNames,ia,ib] = intersect(geneNamesGWAS,geneNamesDrug,'stable');
    geneWeightsGWAS = geneWeightsGWAS(ia);
    drugScores = drugScoresAll(ib,:);

    fprintf(1,'%u matching (/%u %s GWAS); (/%u with treatment weights)\n',...
        length(ia),length(geneWeightsGWAS),whatDisease,length(drugScores));

    %-------------------------------------------------------------------------------
    % Get scores for the property of interest:
    rhos = zeros(numDiseases_Treatment,1);
    for k = 1:numDiseases_Treatment
        geneWeights_treatment = drugScores(:,k);
        rhos(k) = ComputeDotProduct(geneWeights_treatment,geneWeightsGWAS);
        if isnan(rhos(k))
            warning('Issue with %s-%s',whatDiseases_Treatment{k},whatDisease)
        end
    end

    % Generate null distributions:
    numNulls = 1000;
    nullScores = zeros(numNulls,1);
    whatNull = 'randomWeight'; % randomWeight, randomDisease
    for k = 1:numNulls
        switch whatNull
        case 'randomWeight'
            geneWeightsRand = rand(numDrugScores,1);
            geneWeightsRand = geneWeightsRand/norm(geneWeightsRand,2);
            nullScores(k) = ComputeDotProduct(geneWeightsRand,geneWeightsGWAS);
        case 'randomDisease'
            % Shuffle weights taken from a random disease (pooled nulls):
            % [could be done individually for each particular weighting if needed]
            diseaseInd = randi(numDiseases_Treatment,1);
            geneWeightsRand = drugScores(diseaseInd,:)';
            nullScores(k) = ComputeDotProduct(geneWeightsRand,geneWeightsGWAS,true);
        end
    end

    % Compute p-values:
    isSig = zeros(numDiseases_Treatment,1);
    for k = 1:numDiseases_Treatment
        isSig(k) = (mean(rhos(k) < nullScores) < 0.05);
    end

    % Sort:
    [rhos,ix] = sort(rhos,'descend');

    %---------------------------------------------------------------------------
    ax{i} = subplot(1,numDiseases_GWAS,i); hold on
    b = bar(rhos);
    if addNull && ~all(isnan(nullScores))
        % Add null distribution:
        [ff,x] = ksdensity(nullScores,linspace(min(nullScores),max(nullScores),500),'function','pdf');
        ff = 0.8*ff/max(ff);
        plot(ones(2,1)*(numDiseases_Treatment+1)+ff,x,'k');
        plot(ones(2,1)*(numDiseases_Treatment+1)-ff,x,'k');
        ax{i}.XTick = 1:numDiseases_Treatment+1;
        ax{i}.XTickLabel = {whatDiseases_Treatment{ix},'null'};
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
            ax{i}.YLim = [minLim,maxLim];
        end
    else
        ax{i}.XTick = 1:numDiseases_Treatment;
        ax{i}.XTickLabel = whatDiseases_Treatment(ix);
    end
    ax{i}.XTickLabelRotation = 45;
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
linkaxes([ax{:}],'y');
f.Position = [471,945,1830,318];

end
