function DistinguishingCharBar(similarityType,whatProperty, whatNull, whatThreshold)

if nargin < 1
    similarityType = 'MAGMAdefault';
    % {'MAGMAdefault';'Adult_brain';'Fetal_brain';'Neuro';'Astro';'eQTLbrain';'eQTLWhole_Blood';'eQTLLiver';'eQTLHeart_Left_Ventricle';'PPI_mapped_th0';'PPI_eQTLbrain_th0';'PPI_mapped_th400';'PPI_eQTLbrain_th400';'PPI_mapped_th600';'PPI_eQTLbrain_th600';'PPI_mapped_th900';'PPI_eQTLbrain_th900';'AllenMeanCoexpMapped';'AllenMeanCoexpeQTLbrain'}
end
if nargin < 2
    whatProperty = 'P';
    % for MAGMA-based: {'ZSTAT';'P';'NSNPS';'NSNPSnorm'}
    % for PPI-based: {'numPPIneighbors1';'percPPIneighbors1';'weiPPIneighbors1';'expWeiPPIneighbors1';'numPPIneighbors2';'percPPIneighbors2';'weiPPIneighbors2';'expWeiPPIneighbors2';'numPPIneighbors3';'percPPIneighbors3';'weiPPIneighbors3';'expWeiPPIneighbors3';'numPPIneighbors4';'percPPIneighbors4';'weiPPIneighbors4';'expWeiPPIneighbors4';'numPPIneighbors5';'percPPIneighbors5';'weiPPIneighbors5';'expWeiPPIneighbors5';'numPPIneighbors6';'percPPIneighbors6';'weiPPIneighbors6';'expWeiPPIneighbors6';'medianPPIDistance';'meanPPIDistance'}
end
if nargin < 3
    whatNull = 'randomDisease'; 
end
if nargin <4
    whatThreshold = 'BF';
end

addNull = true;

switch whatNull
    case 'randomGene'
        % if using randomGene null, then don't include AD - it has no drug
        % targets
        whatDiseases_GWAS = {'ADHD', 'MDD2', 'SCZ', 'BIP2', 'DIABETES', 'HF'};
    otherwise
        whatDiseases_GWAS = {'ADHD', 'MDD2', 'SCZ', 'BIP2', 'DIABETES', 'HF', 'AD'};
end

whatDiseases_Treatment = {'ADHD','BIP','SCZ','MDD','pulmonary','cardiology','gastro','diabetes'};

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
    [geneNamesGWAS,geneWeightsGWAS] = GiveMeNormalizedScoreVector(whatDisease,'GWAS',similarityType,whatProperty, whatThreshold);

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
    
    if strcmp(whatNull, 'randomWeight') || strcmp(whatNull, 'randomDisease') || strcmp(whatNull, 'randomPsychDisease')
        % one single set of nulls for the whole analysis
        for k = 1:numNulls
            switch whatNull
                case 'randomWeight'
                    geneWeightsRand = rand(numDrugScores,1);
                    nullScores(k) = ComputeDotProduct(geneWeightsRand,geneWeightsGWAS);
                case 'randomDisease'
                    % Shuffle weights taken from a random disease (pooled nulls):
                    % [could be done individually for each particular weighting if needed]
                    diseaseInd = randi(numDiseases_Treatment,1);
                    geneWeightsRand = drugScores(:,diseaseInd);
                    nullScores(k) = ComputeDotProduct(geneWeightsRand,geneWeightsGWAS,true);
                case 'randomPsychDisease'
                    % based on GWAS list, constrain the null to only shuffle the
                    % psychiatric disorders or non-psychiatric;
                    
                    psychDIS = contains(whatDiseases_Treatment, 'ADHD') | contains(whatDiseases_Treatment, 'BIP') | ...
                        contains(whatDiseases_Treatment, 'SCZ') | contains(whatDiseases_Treatment, 'MDD');
                    
                    % find columns for psychiatric drug lists
                    if strcmp(whatDisease, 'ADHD') || strcmp(whatDisease,'MDD2') || ...
                            strcmp(whatDisease, 'SCZ') || strcmp(whatDisease, 'BIP2')
                        selectDIS_IND = find(psychDIS);
                    else
                        selectDIS_IND = find(psychDIS==0);
                    end
                    
                    num_DIS = length(selectDIS_IND);
                    diseaseInd = selectDIS_IND(randi(num_DIS,1));
                    geneWeightsRand = drugScores(:,diseaseInd);
                    nullScores(k) = ComputeDotProduct(geneWeightsRand,geneWeightsGWAS,true);
            end
        end
    else
        
        % separate set of nulls for each drug target list
        switch whatNull
            
            case 'randomGene'
                % select a random set of genes from GWAS scores - keep
                % drugs the same, randomise GWAS scores
                
                nullScores(k,l) = ComputeDotProduct(geneWeightsRand,geneWeightsGWAS_rand,true);
            case 'randomDrug'
                % shuffle genes within the drugs each disorder separately
                % - do drugs for the disorder match GWAS hits better than a
                % random set of gene targets with the same properties;
                % this can only be done for GWAS lists with drugs; , so not for AD;
                switch whatDisease
                    case 'HF'
                        diseaseInd = contains(whatDiseases_Treatment, 'cardiology', 'IgnoreCase',true);
                    otherwise
                        indLET = isletter(whatDisease);
                        diseaseInd = contains(whatDiseases_Treatment, whatDisease(indLET), 'IgnoreCase',true);
                end
                
                drugScores_DIS = drugScores(:,diseaseInd);
                geneWeightsRand = drugScores_DIS(randperm(numel(drugScores_DIS)));
                nullScores(k,l) = ComputeDotProduct(geneWeightsRand,geneWeightsGWAS,true);
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
            if minLim > null_p10
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
