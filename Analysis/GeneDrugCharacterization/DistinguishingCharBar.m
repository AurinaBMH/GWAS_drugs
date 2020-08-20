function [rhosALL ,pValsALL, whatDiseases_Treatment, geneWeights_treatment, geneWeightsGWAS] = DistinguishingCharBar(similarityType,whatProperty, whatNull, whatThreshold, whatDiseases_GWAS, doPlot)

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

if nargin<5
    whatDiseases_GWAS = {'ADHD', 'MDD2', 'SCZ', 'BIP2', 'DIABETES', 'HF', 'AD'};
end

if nargin<6
    doPlot = true;
end

whatDiseases_Treatment = {'ADHD','BIP','SCZ','MDD','pulmonary','cardiology','gastro','diabetes'};
params = SetDefaultParams();
whatScore = params.whatScore;

if strcmp(whatNull, 'randomGene')
    load('GWAS_disordersMAGMA.mat', 'DISORDERlist')
elseif strcmp(whatNull, 'randomDrug')
    load(sprintf('nulls_5000_%stargets_randomDrug.mat', params.whatDrugTargets), 'RANDOMdrugs_treatment', 'whatDiseases_Treatment', 'geneNames');
    geneNames_nulls = geneNames; 
end
addNull = true;

%-------------------------------------------------------------------------------
% Load in default parameters:


%-------------------------------------------------------------------------------
numDiseases_Treatment = length(whatDiseases_Treatment);
numDiseases_GWAS = length(whatDiseases_GWAS);

%-------------------------------------------------------------------------------
% Load treatment weights of each gene implicated in each disorder:
[geneNamesDrug,drugScoresAll] = GiveMeNormalizedScoreVectors(whatDiseases_Treatment,'Drug');
numDrugScores = length(drugScoresAll);

%===============================================================================
if doPlot
    f = figure('color','w');
    sgtitle(sprintf('%s, %s',similarityType, whatProperty))
    ax = cell(numDiseases_GWAS,1);
    if length(whatDiseases_GWAS)>3
        f.Position = [471,945,1830,318];
    end
end
rhosALL = zeros(numDiseases_Treatment,numDiseases_GWAS);
pValsALL = zeros(numDiseases_Treatment,numDiseases_GWAS);

for i = 1:numDiseases_GWAS
    whatDisease = whatDiseases_GWAS{i};
    [geneNamesGWAS,geneWeightsGWAS] = GiveMeNormalizedScoreVector(whatDisease,'GWAS',similarityType,whatProperty, whatThreshold);
    
    % Combine two datasets on overlap:
    [geneNames,ia,ib] = intersect(geneNamesGWAS,geneNamesDrug);
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
    
    rhosALL(:,i) = rhos;
    
    % Generate null distributions:
    numNulls = 1000; 
    isSig = zeros(numDiseases_Treatment,1);
    pVals = zeros(numDiseases_Treatment,1);
    all_nullScores = cell(numDiseases_Treatment,1);
    
    if strcmp(whatNull, 'randomWeight') || strcmp(whatNull, 'randomDisease') || strcmp(whatNull, 'randomPsychDisease')
        % one single set of nulls for the whole analysis
        nullScores = zeros(numNulls,1);
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
                    % Shuffle weights taken from a random
                    % psychiatric disorders or non-psychiatric disease (pooled nulls):
                    
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
        % Compute p-values based on pooled nulls
        for k = 1:numDiseases_Treatment
            isSig(k) = (mean(rhos(k) < nullScores) < 0.05);
            pVals(k) = mean(rhos(k) < nullScores);
        end
        
    else
        
        for l = 1:numDiseases_Treatment
            nullScores = zeros(numNulls,1);
            for k = 1:numNulls
                % separate set of nulls for each drug target list
                switch whatNull
                    
                    case 'randomGene' % is the actual match higher than a match with completely random genes
                        if ~contains(similarityType, 'PPI') && ~contains(similarityType, 'Allen')
                            % for this null, load all available scores for genes
                            % select a random set of genes from GWAS scores - keep
                            % drugs the same, randomise GWAS scores; This is
                            % suitable only for MAGMA-based  methods;
                            
                            switch whatProperty
                                case 'P'
                                    geneWeightsGWAS_all = -log10(DISORDERlist.(similarityType).(whatDisease).(whatProperty));
                                otherwise
                                    geneWeightsGWAS_all = DISORDERlist.(similarityType).(whatDisease).(whatProperty);
                            end
                            
                            geneWeightsGWAS_rand = datasample(geneWeightsGWAS_all,numDrugScores,'Replace',false);
                            % normalizde the weights, by default this used
                            % norm-1, but mabe should be chnaged to norm2?
                            geneWeightsGWAS_randNorm = normalizeScoreVector(geneWeightsGWAS_rand);
                            nullScores(k) = ComputeDotProduct(drugScores(:,l),geneWeightsGWAS_randNorm);
                        else
                            
                            warning('% null is not compatible with %s\n', whatNull, whatNull)
                            
                        end
                        
                    case 'randomTarget' % is the actual match higher than a match with random gene score assignment
                        % Shuffle weights taken from each drug list individually
                        drugScores_DIS = drugScores(:,l);
                        nullScores(k) = ComputeDotProduct(drugScores_DIS,geneWeightsGWAS, true);
                        % randomise v1 within ComputeDotProduct
                    case 'randomDrug' % for each disease get a random set of drugs that is the same size as
                        % real list of drugs, e.g. for ADHD select 18 drugs
                        % and get a normalized score vector for this list of drugs as if it's a separate disease
                        %drugScores_DIS = give_randomDrug_null(whatDiseases_Treatment{l}, disorderDrugs, allDrugs);
                        % load pre-computed nulls
                        [~, ix, iy] = intersect(geneNames, geneNames_nulls); 
                        drugScores_DIS = RANDOMdrugs_treatment{l}(iy,k);
                        %drugScores_DIS = RANDOMdrugs_treatment{l}(:,k);
                        nullScores(k) = ComputeDotProduct(drugScores_DIS,geneWeightsGWAS(ix));
                        %nullScores(k) = ComputeDotProduct(drugScores_DIS,geneWeightsGWAS(:));
                        
                end
            end
            % Compute p-values: based on separate nulls
            isSig(l) = (mean(rhos(l) < nullScores) < 0.05);
            pVals(l) = mean(rhos(l) < nullScores);
            % save separate nulls
            all_nullScores{l} = nullScores;
        end
    end
    
    pValsALL(:,i) = pVals;
    
    % Sort:
    [rhos,ix] = sort(rhos,'descend');
    all_nullScores = all_nullScores(ix); 
    
    %---------------------------------------------------------------------------
    if doPlot
        ax{i} = subplot(1,numDiseases_GWAS,i); hold on
        b = bar(rhos);
        
        if strcmp(whatNull, 'randomWeight') || strcmp(whatNull, 'randomDisease') || strcmp(whatNull, 'randomPsychDisease')
            if addNull && ~all(isnan(nullScores))
                
                % function to plot nulls
                % this will plot the distribution at the end of bar chart for pulled nulls
                [minLim,maxLim] = plot_nullDistribution(nullScores, rhos, numDiseases_Treatment+1);
                
                ax{i}.XTick = 1:numDiseases_Treatment+1;
                ax{i}.XTickLabel = {whatDiseases_Treatment{ix},'null'};
                
                if range(rhos) > 0
                    ax{i}.YLim = [minLim,maxLim];
                end
                
            else
                
                ax{i}.XTick = 1:numDiseases_Treatment;
                ax{i}.XTickLabel = whatDiseases_Treatment(ix);
                
            end
            
        else
            
            % for other nulls plot distributions on each bar
            minLim = zeros(numDiseases_Treatment,1); 
            maxLim = zeros(numDiseases_Treatment,1); 
            
            for nn = 1:numDiseases_Treatment
                nullScores_selected = all_nullScores{nn}; 
                if addNull && ~all(isnan(nullScores_selected))
                    
                    % function to plot nulls
                    % this will plot the distribution at the end of baf chart
                    % all_nullScores are now reordered to correspond ro rhos 
                    [minLim(nn),maxLim(nn)] = plot_nullDistribution(nullScores_selected, rhos, nn);
                    
                end
 
            end
            
            if range(rhos) > 0
                ax{i}.YLim = [min(minLim),max(maxLim)];
            end
            ax{i}.XTick = 1:numDiseases_Treatment;
            ax{i}.XTickLabel = whatDiseases_Treatment(ix);
            
        end
        
        
        ax{i}.XTickLabelRotation = 45;
        xlabel('Disease treatment')
        ylabel(sprintf('%s similarity',whatScore))
        title(sprintf('%s-%s',whatDiseases_GWAS{i},whatProperty),'interpreter','none')
        cMapGeneric = BF_getcmap('set2',numDiseases_Treatment,false);
        for k = find(~isSig)'
            cMapGeneric(k,:) = brighten(cMapGeneric(k,:),0.95);
        end
        b.CData = cMapGeneric(ix,:);
        b.FaceColor = 'flat';
        
        linkaxes([ax{:}],'y');
    end
    
    
end
end
