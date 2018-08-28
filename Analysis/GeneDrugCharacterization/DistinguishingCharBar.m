function DistinguishingCharBar(whatProperty)

if nargin < 1
    whatProperty = 'numGWASMapped';
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
% indicatorTable = percIndicatorTable;

%===============================================================================
f = figure('color','w');
for i = 1:numDiseases_GWAS
    % Load data:
    fileName = sprintf('resultsTable_%s.mat',whatDiseases_GWAS{i});
    load(fileName,'geneScores');
    switch whatProperty
    case 'numGWASMapped'
        geneWeights_GWAS = geneScores.DNA.numGWAS;
    case 'percGWASMapped'
        geneWeights_GWAS = geneScores.DNA.percGWAS;
    case 'PPIMappedWeighted_percNeigh1'
        geneWeights_GWAS = geneScores.PPI_mapped_weighted.percPPIneighbors1;
    case 'PPI_mapped_th0_percNeigh1'
        geneWeights_GWAS = geneScores.PPI_mapped_th0.percPPIneighbors1;
    case 'PPI_mapped_th400_percNeigh1'
        geneWeights_GWAS = geneScores.PPI_mapped_th400.percPPIneighbors1;
    case 'PPI_mapped_th400_percNeigh2'
        geneWeights_GWAS = geneScores.PPI_mapped_th400.percPPIneighbors2;
    case 'PPI_mapped_th400_numNeigh2'
        geneWeights_GWAS = geneScores.PPI_mapped_th400.numPPIneighbors2;
    case 'PPI_mapped_th400_percNeigh3'
        geneWeights_GWAS = geneScores.PPI_mapped_th400.percPPIneighbors3;
    case 'PPI_mapped_th400_percNeigh4'
        geneWeights_GWAS = geneScores.PPI_mapped_th400.percPPIneighbors4;
    otherwise
        error('Unknown property: ''%s''',whatProperty);
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
            geneWeights_treatment = indicatorTable.(whatDiseases_Treatment{k});
            % We want the treatment vector to sum to 1:
            geneWeights_treatment = geneWeights_treatment/sum(geneWeights_treatment);
            r = ~isnan(geneWeights_treatment) & ~isnan(geneWeights_GWAS);
            rhos(k) = sum(geneWeights_treatment(r).*geneWeights_GWAS(r));
        end

        % Generate (pooled) nulls [could be done individually for each particular weighting if needed]:
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
        if min(rhos) < max(rhos)
            ax.YLim = [min(rhos)*0.9,max(rhos)*1.1];
        end
    else
        ax.XTick = 1:numDiseases_Treatment;
        ax.XTickLabel = whatDiseases_Treatment(ix);
    end
    ax.XTickLabelRotation = 45;
    xlabel('Disease treatment')
    ylabel(whatScore)
    title(sprintf('%s-%s',whatDiseases_GWAS{i},whatProperty),'interpreter','none')
    cMapGeneric = BF_getcmap('set2',numDiseases_Treatment,false);
    for k = find(~isSig)'
        cMapGeneric(k,:) = brighten(cMapGeneric(k,:),0.8);
    end
    b.CData = cMapGeneric(ix,:);
    b.FaceColor = 'flat';
    %-------------------------------------------------------------------------------
    % Add null distribution

end
