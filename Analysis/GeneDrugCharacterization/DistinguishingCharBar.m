function [rhosALL ,pValsALL, whatDiseases_Treatment, geneWeights_treatment, geneWeightsGWAS, nullScoresALL] = DistinguishingCharBar(similarityType,whatProperty, whatNull, whatThreshold, whatDiseases_GWAS, doPlot, numMeasures, whatMeasures)

if nargin < 1
    similarityType = 'MAGMAdefault';
end
if nargin < 2
    whatProperty = 'P';
end
if nargin < 3
    whatNull = 'randomDrugR_all_drugbank';
end
if nargin < 4
    whatThreshold = 'BF';
end

if nargin < 5
    params = SetDefaultParams();
    whatDiseases_GWAS = params.whatGWAS;
end

if nargin < 6
    doPlot = true;
end

if nargin < 7
    numMeasures = length(similarityType); 
end

%-------------------------------------------------------------------------------
% Load in default parameters:
params = SetDefaultParams();    
whatDiseases_Treatment_ALL = params.whatDiseases_Treatment_ALL;
whatScore = params.whatScore;
%-------------------------------------------------------------------------------

switch whatMeasures
    case {'allPsych', 'reduced'}

        whatDiseases_Treatment_SEL = params.whatDiseases_Treatment; 
        whatDiseases_Treatment_label = params.whatDiseases_Treatment_label; 
    case 'allBody'

        whatDiseases_Treatment_SEL = params.whatDiseases_Treatment_body; 
        whatDiseases_Treatment_label = params.whatDiseases_Treatment_label_body;
        
    case 'all'

        whatDiseases_Treatment_SEL = params.whatDiseases_Treatment_label_ALL; 
        whatDiseases_Treatment_label = params.whatDiseases_Treatment_label_ALL; 
        
end

if strcmp(whatNull, 'randomGene')
    
    load('GWAS_disordersMAGMA.mat', 'DISORDERlist')
    
elseif strcmp(whatNull, 'randomDrugR') || strcmp(whatNull, 'randomDrugP') || ...
        strcmp(whatNull, 'randomDrugP_all_drugbank_psych') || strcmp(whatNull, 'randomDrugR_all_drugbank') || ...
        strcmp(whatNull, 'randomDrugP_active_drugbank_psych') || strcmp(whatNull, 'randomDrugR_active_drugbank')
        
    
    load(sprintf('nulls_5000_%stargets_%s.mat', params.whatDrugTargets, whatNull), 'RANDOMdrugs_treatment', 'whatDiseases_Treatment', 'geneNames');
    %geneNames_nulls = geneNames;
    % select nulls for drugs that will be visualised
    [whatDiseases_Treatment, Tind] = intersect(whatDiseases_Treatment, whatDiseases_Treatment_SEL, 'stable');
    RANDOMdrugs_treatment = RANDOMdrugs_treatment(Tind);
     
else
    
    whatDiseases_Treatment = whatDiseases_Treatment_SEL;
    
end
% select only relevant nulls
addNull = true;

numDiseases_Treatment = length(whatDiseases_Treatment);
numDiseases_GWAS = length(whatDiseases_GWAS);

%-------------------------------------------------------------------------------
% Load treatment weights of each gene implicated in each disorder:
% load drugs scores for all disorders, will select values from these
[geneNamesDrug,drugScoresAll] = GiveMeNormalizedScoreVectors(whatDiseases_Treatment_ALL,'Drug');
numDrugScores = length(drugScoresAll);

%===============================================================================
if doPlot
    f = figure('color','w', 'Position', [300, 300, 1500, 400]);
    %sgtitle(sprintf('%s, %s',similarityType, whatProperty))
    ax = cell(numDiseases_GWAS,1);
    if length(whatDiseases_GWAS)>3
        f.Position = [471,945,1830,318];
    end
end
rhosALL = zeros(numDiseases_Treatment,numDiseases_GWAS);
pValsALL = zeros(numDiseases_Treatment,numDiseases_GWAS);
nullScoresALL = cell(numDiseases_GWAS,1); 
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
        % select drug list
        kIND = contains(whatDiseases_Treatment_ALL, whatDiseases_Treatment{k});
        geneWeights_treatment = drugScores(:,kIND);
        
        rhos(k) = ComputeDotProduct(geneWeights_treatment,geneWeightsGWAS);
        if isnan(rhos(k))
            warning('Issue with %s-%s',whatDiseases_Treatment{k},whatDisease)
        end
    end
    
    rhosALL(:,i) = rhos;
    
    % Generate null distributions:
    numNulls = params.numNull;
    isSig_BF = false(numDiseases_Treatment,1);
    isSig = false(numDiseases_Treatment,1);
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
            isSig_BF(k) = logical((mean(rhos(k) < nullScores) < 0.05/numMeasures));
            isSig(k) = logical((mean(rhos(k) < nullScores) < 0.05));
            pVals(k) = mean(rhos(k) < nullScores);
        end
        
    else
        
        for l = 1:numDiseases_Treatment
            % select drugs from selected list
            lIND = contains(whatDiseases_Treatment_ALL, whatDiseases_Treatment{l});
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
                            nullScores(k) = ComputeDotProduct(drugScores(:,lIND),geneWeightsGWAS_randNorm);
                        else
                            
                            warning('% null is not compatible with %s\n', whatNull, whatNull)
                            return
                            
                        end
                        
                    case 'randomTarget' % is the actual match higher than a match with random gene score assignment
                        % Shuffle weights taken from each drug list individually
                        drugScores_DIS = drugScores(:,lIND);
                        nullScores(k) = ComputeDotProduct(drugScores_DIS,geneWeightsGWAS, true);
                        % randomise v1 within ComputeDotProduct
                    case {'randomDrugP','randomDrugR',...
                            'randomDrugP_all_drugbank_psych', 'randomDrugR_all_drugbank', ...
                            'randomDrugP_active_drugbank_psych', 'randomDrugR_active_drugbank'}  % for each disease get a random set of drugs that is the same size as
                        % real list of drugs, e.g. for ADHD select 18 drugs
                        % load pre-computed nulls
                        % the order of genes is the same in GWAS and nulls
                        % [~, ix, iy] = intersect(geneNames, geneNames_nulls);
                        %  drugScores_DIS = RANDOMdrugs_treatment{l}(iy,k);
                        
                        % using l, because RANDOMdrugs_treatment is already
                        % reduced only to drugs that are being tested
                        drugScores_DIS = RANDOMdrugs_treatment{l}(:,k);
                        %nullScores(k) = ComputeDotProduct(drugScores_DIS,geneWeightsGWAS(ix));
                        nullScores(k) = ComputeDotProduct(drugScores_DIS,geneWeightsGWAS);
                        
                end
            end
            % Compute p-values: based on separate nulls
            isSig_BF(l) = logical((mean(rhos(l) < nullScores) < 0.05/numMeasures));
            isSig(l) = logical((mean(rhos(l) < nullScores) < 0.05));
            pVals(l) = mean(rhos(l) < nullScores);
            % save separate nulls
            all_nullScores{l} = nullScores;
            
        end
    end
    
    pValsALL(:,i) = pVals;
    
    
    % Sort:
    %[~, ix] = sort(pVals, 'ascend');
    [~,ix] = sort(rhos,'descend');
    nullScoresALL{i} = all_nullScores; 
    all_nullScores = all_nullScores(ix);
    rhos = rhos(ix);
    
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
                ax{i}.XTickLabel = {whatDiseases_Treatment_label{ix},'null'};
                
                if range(rhos) > 0
                    ax{i}.YLim = [minLim,maxLim];
                end
                
            else
                
                ax{i}.XTick = 1:numDiseases_Treatment;
                ax{i}.XTickLabel = whatDiseases_Treatment_label(ix);
                
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
            ax{i}.XTickLabel = whatDiseases_Treatment_label(ix);
            
        end
        
        
        ax{i}.XTickLabelRotation = 45;
        xlabel('Disease treatment')
        %ylabel(sprintf('%s similarity',whatScore))
        ylabel({'GWAS-treatment', 'similarity'})
        %title(sprintf('%s-%s',whatDiseases_GWAS{i},whatProperty),'interpreter','none')
        title(sprintf('%s',whatDiseases_GWAS{i}),'interpreter','none')
        cMapGeneric = BF_getcmap('set4',numDiseases_Treatment,false);
        cMapGeneric_n = cMapGeneric;
        
        for k=1:length(isSig)
            if ~isSig(k)
                cMapGeneric_n(k,:) = brighten(cMapGeneric(k,:),0.99);
            elseif ~isSig_BF(k) && isSig(k)
                cMapGeneric_n(k,:) = brighten(cMapGeneric(k,:),0.85);
            end
        end

        b.CData = cMapGeneric_n(ix,:);
        b.FaceColor = 'flat';
        
        linkaxes([ax{:}],'y');
        set(gca,'FontSize', 14)
    end
    
    
end
end
