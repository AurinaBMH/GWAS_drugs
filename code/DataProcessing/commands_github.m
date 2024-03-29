% binary full: PPI Adj and distance matrix
[AdjPPI,geneNames] = PPINImport(false,0,'HGNCmatch');
distMatrix = ComputePPIDist(0,false);

% binary >400 evidence score: PPI Adj and distance matrix
[AdjPPI,geneNames] = PPINImport(false,400,'HGNCmatch');
distMatrix = ComputePPIDist(400,false);

% binary >600 evidence score: PPI Adj and distance matrix
[AdjPPI,geneNames] = PPINImport(false,600,'HGNCmatch');
distMatrix = ComputePPIDist(600,false);

% binary >900 evidence score: PPI Adj and distance matrix
[AdjPPI,geneNames] = PPINImport(false,900,'HGNCmatch');
distMatrix = ComputePPIDist(900,false);

% weighted
[AdjPPI,geneNames] = PPINImport(true);


distMatrix = ComputePPIDist([],true); - % this takes ages, didn't finish in 3 days.

% now result tables are generated for each disorder
GenerateResultsTables;

% loop over all mapping methods and save the results to file, so can compare;
% load ADHD as example to select mapping methods
load('resultsTable_ADHD_BF_2020.mat', 'geneScores')
similarityTypes = setdiff(fieldnames(geneScores), {'gene', 'params',...
    'PPI_eQTLbrain_th0', 'PPI_eQTLbrain_th400',...
    'PPI_mapped_th0', 'PPI_mapped_th400',...
    'AllenMeanCoexpMapped', 'AllenMeanCoexpeQTLbrain'});

whatThreshold = 'BF';
whatNull = 'randomDrug';

for t=1:length(similarityTypes)

    if contains(similarityTypes{t},'PPI')

        whatProperty = 'percPPIneighbors1';

    elseif contains(similarityTypes{t},'Allen')
        whatProperty = 'r';
    else
        whatProperty = 'P';
    end

    DistinguishingCharBar(similarityTypes{t},whatProperty, whatNull, whatThreshold);
    figureName = sprintf('figures/GWASdrug_%s_%s_%s_%s', similarityTypes{t}, whatProperty, whatNull, whatThreshold);
    print(gcf,figureName,'-dpng','-r300');


end

% visualise GWAS vs drug targetsas scatter
whatDiseaseGWAS = 'ADHD';
whatProperty = 'MAGMAdefault';
whatDiseaseDrug = 'ADHD';
whatThreshold = 'BF';

GiveMeComparisonScores(whatDiseaseGWAS,whatProperty,whatDiseaseDrug,whatThreshold)



% visualise gene score vectors side-by-side for GWAS targets
whichDiseases = {'ADHD','MDD2','BIP2','SCZ','DIABETES','HF', 'AD'};
whatMeasurement = 'GWAS';
similarityTypes = {'MAGMAdefault','eQTLbrain', 'PPI_mapped_th600'};
whatPropertys = {'P', 'P', 'percPPIneighbors1'};

for k=1:length(similarityTypes)

    VisualizeWeightingVectors(whichDiseases, whatMeasurement, similarityTypes{k}, whatPropertys{k});
    figureName = sprintf('figures/geneWeights_%s_%s_%s', whatMeasurement, similarityTypes{k}, whatPropertys{k});
    print(gcf,figureName,'-dpng','-r300');

end


% visualise gene score vectors side-by-side for DRUG lists
whichDiseases = {'ADHD','BIP','SCZ','MDD','pulmonary','cardiology','gastro','diabetes'};
whatMeasurement = 'Drug';
VisualizeWeightingVectors(whichDiseases, whatMeasurement)
figureName = sprintf('figures/geneWeights_%s', whatMeasurement);
print(gcf,figureName,'-dpng','-r300');

% plot score distributions for different disorders
figure;
subplot(1,4,1); histogram(indicatorTable.ADHD); title('ADHD'); xlim([0 0.25]); ylim([0 600])
subplot(1,4,2); histogram(indicatorTable.diabetes); title('DIABETES'); xlim([0 0.25]); ylim([0 600])
subplot(1,4,3); histogram(indicatorTable.cardiology); title('CARDIOLOGY'); xlim([0 0.25]); ylim([0 600])
subplot(1,4,4); histogram(indicatorTable.gastro); title('GASTRO'); xlim([0 0.25]); ylim([0 600])

% visualise GWAS and drug vectors side by side for a pair of selected disorders
% what pairs to visualise:
whatDiseaseGWAS = {'ADHD', 'MDD2', 'SCZ', 'BIP2', 'DIABETES', 'HF'}; % select from: 'ADHD', 'MDD2', 'SCZ', 'BIP2', 'DIABETES', 'HF', 'AD';
whatDiseaseDRUG = {'ADHD', 'MDD', 'SCZ', 'BIP', 'diabetes', 'cardiology'}; % select from: 'ADHD','BIP','SCZ','MDD','pulmonary','cardiology','gastro','diabetes'
similarityType = {'MAGMAdefault', 'PPI_mapped_th600'};
whatProperty = {'P', 'percPPIneighbors1'};


for D = 1:length(whatDiseaseGWAS)
figure('color','w');
    for T = 1:length(similarityType)
        subplot(length(similarityType),1,T);
        data = VisualizeWeightingVectors_GWASvsDRUG(whatDiseaseGWAS{D}, whatDiseaseDRUG{D}, similarityType{T}, whatProperty{T}, false);
        figureName = sprintf('figures/geneWeights_GWASvsDRUG_%s_%s', whatDiseaseGWAS{D}, whatDiseaseDRUG{D});
        print(gcf,figureName,'-dpng','-r300');
    end
end
