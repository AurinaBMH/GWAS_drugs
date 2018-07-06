whatDiseases_GWAS = {'ADHD','BIP','SZP','MDD'};
whatProperty = 'percPPIneighbors1DiseaseLD';
whatDiseases_Treatment = {'ADHD','BIP','SZP','MDD','pulmonary','cardiology','gastro','diabetes'};

%-------------------------------------------------------------------------------
numDiseases_Treatment = length(whatDiseases_Treatment);
numDiseases_GWAS = length(whatDiseases_GWAS);

%===============================================================================
f = figure('color','w');
for i = 1:numDiseases_GWAS
    %-------------------------------------------------------------------------------
    % Load data:
    %-------------------------------------------------------------------------------
    % (1) Similarity of each gene to GWAS hits for this disorder:
    fileName = sprintf('resultsTable_%s-0.0.mat',whatDiseases_GWAS{i});
    load(fileName,'resultsTable');

    % (2) Treatment weights of each gene implicated in the disorder:
    normalizeWithinDrugs = true; % weight genes lower if they occur in drugs with large numbers of gene targets
    [indicatorTable,percIndicatorTable] = ImportTreatmentLists(normalizeWithinDrugs);
    % indicatorTable = percIndicatorTable;

    % Combine two datasets on overlap:
    [~,ia,ib] = intersect(resultsTable.gene,indicatorTable.Properties.RowNames,'stable');
    indicatorTable = indicatorTable(ib,:);
    resultsTable = resultsTable(ia,:);

    %-------------------------------------------------------------------------------
    % Get scores for the property of interest:
    rhos = zeros(numDiseases,1);
    isSig = zeros(numDiseases,1);
    geneWeights_GWAS = resultsTable.(whatProperty);
    for k = 1:numDiseases
        geneWeights_treatment = indicatorTable.(whatDiseases{k});
        [rhos(k),p] = corr(geneWeights_GWAS,geneWeights_treatment,'type','Kendall','rows','pairwise');
        isSig(k) = (p < 0.05);
    end
    [rhos,ix] = sort(rhos,'descend');

    subplot(2,2,i)
    b = bar(rhos);
    ax = gca;
    ax.XTick = 1:numDiseases;
    ax.XTickLabel = whatDiseases(ix);
    xlabel('Disease treatment')
    ylabel('Kendell-corr')
    title(sprintf('%s-%s',whatDiseases_GWAS{i},whatProperty))
    cMapGeneric = BF_getcmap('set2',numDiseases,false);
    for k = find(~isSig)'
        cMapGeneric(k,:) = brighten(cMapGeneric(k,:),0.8);
    end
    b.CData = cMapGeneric(ix,:);
    b.FaceColor = 'flat';

end
