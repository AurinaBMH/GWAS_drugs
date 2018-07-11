function DistinguishingCharBar(whatProperty)

if nargin < 1
    whatProperty = 'numGWASMapped';
    % 'numGWASMapped','numLDSNPs','percPPIneigh1Mapped','percPPIneigh1LD','AllenMeanCoexpMapped','AllenMeanCoexpLD'
end

whatScore = 'weightedSum'; %'Kendall', 'weightedSum'
whatDiseases_GWAS = {'ADHD','BIP','SZP','MDD','diabetes'};
whatDiseases_Treatment = {'ADHD','BIP','SZP','MDD','pulmonary','cardiology','gastro','diabetes'};

%-------------------------------------------------------------------------------
numDiseases_Treatment = length(whatDiseases_Treatment);
numDiseases_GWAS = length(whatDiseases_GWAS);

%-------------------------------------------------------------------------------
% Load treatment weights of each gene implicated in each disorder:
normalizeWithinDrugs = true; % weight genes lower if they occur in drugs with large numbers of gene targets
[indicatorTable,percIndicatorTable] = ImportTreatmentLists(normalizeWithinDrugs);
% indicatorTable = percIndicatorTable;

%===============================================================================
f = figure('color','w');
for i = 1:numDiseases_GWAS
    %-------------------------------------------------------------------------------
    % Load data:
    %-------------------------------------------------------------------------------
    % (1) Similarity of each gene to GWAS hits for this disorder:
    fileName = sprintf('resultsTable_%s-0.0.mat',whatDiseases_GWAS{i});
    load(fileName,'resultsTable');
    geneWeights_GWAS = resultsTable.(whatProperty);

    % Combine two datasets on overlap:
    [~,ia,ib] = intersect(resultsTable.gene,indicatorTable.Properties.RowNames,'stable');
    indicatorTable = indicatorTable(ib,:);
    resultsTable = resultsTable(ia,:);

    fprintf(1,'%u matching (/%u %s GWAS); (/%u with treatment weights)\n',...
            length(ia),height(resultsTable),whatDiseases_GWAS{i},height(indicatorTable));

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
            geneWeights_treatment = indicatorTable.(whatDiseases_Treatment{k});
            % We want the treatment vector to sum to 1:
            geneWeights_treatment = geneWeights_treatment/sum(geneWeights_treatment);
            r = ~isnan(geneWeights_treatment) & ~isnan(geneWeights_GWAS);
            rhos(k) = sum(geneWeights_treatment(r).*geneWeights_GWAS(r));
        end
        % Generate nulls:
        numNulls = 2000;
        nullScores = zeros(numNulls,1);
        for k = 1:numNulls
            rp = randperm(numDiseases_Treatment);
            geneWeights_treatment = indicatorTable.(whatDiseases_Treatment{rp(1)});
            rp = randperm(length(geneWeights_treatment));
            % We want the treatment vector to sum to 1:
            geneWeights_treatment = geneWeights_treatment(rp)/sum(geneWeights_treatment);
            r = ~isnan(geneWeights_treatment) & ~isnan(geneWeights_GWAS);
            nullScores(k) = sum(geneWeights_treatment(r).*geneWeights_GWAS(r));
        end
        % Compute p-values:
        for k = 1:numDiseases_Treatment
            isSig(k) = (mean(rhos(k)<nullScores) < 0.05);
        end
    end
    [rhos,ix] = sort(rhos,'descend');

    %---------------------------------------------------------------------------
    ax = subplot(2,3,i); hold on
    b = bar(rhos);
    if addNull
        [f,x] = ksdensity(nullScores,linspace(min(nullScores),max(nullScores),500),'function','pdf');
        f = f/max(f);
        plot(ones(2,1)*(numDiseases_Treatment+1)+f,x,'k');
        plot(ones(2,1)*(numDiseases_Treatment+1)-f,x,'k');
        ax.XTick = 1:numDiseases_Treatment+1;
        ax.XTickLabel = {whatDiseases_Treatment{ix},'null'};
        ax.YLim = [min(rhos),max(rhos)*1.1];
    else
        ax.XTick = 1:numDiseases_Treatment;
        ax.XTickLabel = whatDiseases_Treatment(ix);
    end
    ax.XTickLabelRotation = 45;
    xlabel('Disease treatment')
    ylabel(whatScore)
    title(sprintf('%s-%s',whatDiseases_GWAS{i},whatProperty))
    cMapGeneric = BF_getcmap('set2',numDiseases_Treatment,false);
    for k = find(~isSig)'
        cMapGeneric(k,:) = brighten(cMapGeneric(k,:),0.8);
    end
    b.CData = cMapGeneric(ix,:);
    b.FaceColor = 'flat';
    %-------------------------------------------------------------------------------
    % Add null distribution

end
