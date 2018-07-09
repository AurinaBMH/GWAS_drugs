function DistinguishingCharBar(whatProperty)

if nargin < 1
    whatProperty = 'percPPIneigh1Mapped';
    % 'numGWASMapped','numLDSNPs','percPPIneigh1Mapped','percPPIneigh1LD','AllenMeanCoexpMapped','AllenMeanCoexpLD'
end

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
    keyboard
    [~,ia,ib] = intersect(resultsTable.gene,indicatorTable.Properties.RowNames,'stable');
    indicatorTable = indicatorTable(ib,:);
    resultsTable = resultsTable(ia,:);

    fprintf(1,'%u matching (/%u %s GWAS); (/%u with treatment weights)\n',...
            length(ia),height(resultsTable),whatDiseases_GWAS{i},height(indicatorTable));

    %-------------------------------------------------------------------------------
    % Get scores for the property of interest:
    rhos = zeros(numDiseases_Treatment,1);
    isSig = zeros(numDiseases_Treatment,1);
    for k = 1:numDiseases_Treatment
        geneWeights_treatment = indicatorTable.(whatDiseases_Treatment{k});
        [rhos(k),p] = corr(geneWeights_GWAS,geneWeights_treatment,...
                                'type','Kendall','rows','pairwise');
        isSig(k) = (p < 0.05);
    end
    [rhos,ix] = sort(rhos,'descend');

    subplot(2,3,i)
    b = bar(rhos);
    ax = gca;
    ax.XTick = 1:numDiseases_Treatment;
    ax.XTickLabel = whatDiseases_Treatment(ix);
    ax.XTickLabelRotation = 45;
    xlabel('Disease treatment')
    ylabel('Kendell-corr')
    title(sprintf('%s-%s',whatDiseases_GWAS{i},whatProperty))
    cMapGeneric = BF_getcmap('set2',numDiseases_Treatment,false);
    for k = find(~isSig)'
        cMapGeneric(k,:) = brighten(cMapGeneric(k,:),0.8);
    end
    b.CData = cMapGeneric(ix,:);
    b.FaceColor = 'flat';

end
