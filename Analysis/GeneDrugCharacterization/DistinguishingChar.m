function DistinguishingChar(whatDisease,resultsTable)
% What features distinguish actual drugs being used?
%-------------------------------------------------------------------------------
if nargin < 1
    whatDisease = 'SZP';
end

% Load data:
if nargin < 2 || isempty(resultsTable)
    % Uses default PPIN threshold (0.0):
    fileName = sprintf('resultsTable_%s-0.4.mat',whatDisease);
    load(fileName,'resultsTable');
    % resultsTable = pipeline(whatDisease);
end

% Need to define a set of properties to compare
propsToCompare = {'numGWASMapped','numLDSNPs','percPPIneighbors1DiseaseMapped',...
            'percPPIneighbors1DiseaseLD','AllenMeanCoexp'};

%-------------------------------------------------------------------------------
% Ok, so now we load information on drugs already used to treat the disorder:
% drugTable = ImportDrugsTreatment();
%
% switch whatDisease
% case 'SZP'
%     fprintf(1,'Using drugs labeled for Schizophrenia\n');
%     theGenesTreat = unique(drugTable.Protein_DrugIDs(drugTable.isSchiz));
%     % fprintf(1,'Using drugs labeled for N05A: Antipsychotics\n');
%     % theGenesTreat = unique(drugTable.Protein_DrugIDs(drugTable.N05A_Antipsychotics));
% case 'ADHD'
%     theGenesTreat = unique(drugTable.Protein_DrugIDs(drugTable.isADHD));
% end
% numGenesTreat = length(theGenesTreat);
%-------------------------------------------------------------------------------
% Now we look for differences in these columns compared to others:
% isTreated = ismember(resultsTable.gene,theGenesTreat);
% fprintf(1,'%u/%u ''treatment genes'' match GWAS-based analysis\n',...
%                             sum(isTreated),numGenesTreat);

%-------------------------------------------------------------------------------
% Get genes for drugs for all diseases of interest:
normalizeWithinDrugs = true; % weight genes lower if they occur in drugs with large numbers of gene targets
[indicatorTable,percIndicatorTable] = ImportTreatmentLists(normalizeWithinDrugs);

whatDiseases = {'ADHD', 'BIP', 'SCZ', 'MDD', 'DIABETES', 'IBD', 'HF', 'RA', 'gastro', 'pulmonary'};
numDiseases = length(whatDiseases);

%-------------------------------------------------------------------------------
% Reorder by resultsTable (so both are matched on gene):
[~,ia,ib] = intersect(resultsTable.gene,indicatorTable.Properties.RowNames,'stable');
indicatorTable = indicatorTable(ib,:);
resultsTable = resultsTable(ia,:);

%-------------------------------------------------------------------------------
% Compare columns in resultsTable
numProps = length(propsToCompare);
f = figure('color','w');
for i = 1:numProps
    ax = subplot(2,3,i);
    dataVector = resultsTable.(propsToCompare{i});
    dataCell = cell(numDiseases,1);
    for k = 1:numDiseases
        geneWeights = indicatorTable.(whatDiseases{k});
        dataCell{k} = dataVector.*geneWeights;
        % IGNORE any NaNs:
        dataCell{k} = dataCell{k}(~isnan(dataCell{k}));
    end
    % Reorder by mean:
    means = cellfun(@mean,dataCell);
    [~,ix] = sort(means,'descend');
    extraParams = struct();
    extraParams.theColors = BF_getcmap('set2',numDiseases,1);
    extraParams.theColors = extraParams.theColors(ix,:);
    [ff,xx] = BF_JitteredParallelScatter(dataCell(ix),true,true,false,extraParams);
    ax.XTick = 1:numDiseases;
    ax.XTickLabels = whatDiseases(ix);
    title(sprintf('%s-%s',whatDisease,propsToCompare{i}));
    xlabel('Drug class')
end

end
