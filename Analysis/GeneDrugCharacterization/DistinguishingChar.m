function DistinguishingChar(whatDisease,resultsTable)
% What features distinguish actual drugs being used?
%-------------------------------------------------------------------------------
if nargin < 1
    whatDisease = 'SZP';
end

% Load data:
if nargin < 2 || isempty(resultsTable)
    % Uses default PPIN threshold:
    resultsTable = pipeline(whatDisease);
end

%-------------------------------------------------------------------------------
% Ok, so now we load information on drugs already used to treat the disorder:
drugTable = ImportDrugsTreatment();

switch whatDisease
case 'SZP'
    fprintf(1,'Using drugs labeled for Schizophrenia\n');
    theGenesTreat = unique(drugTable.Protein_DrugIDs(drugTable.isSchiz));
    % fprintf(1,'Using drugs labeled for N05A: Antipsychotics\n');
    % theGenesTreat = unique(drugTable.Protein_DrugIDs(drugTable.N05A_Antipsychotics));
case 'ADHD'
    theGenesTreat = unique(drugTable.Protein_DrugIDs(drugTable.isADHD));
end
numGenesTreat = length(theGenesTreat);
fprintf(1,'%u genes are targeted by existing drugs for %s:\n',numGenesTreat,whatDisease);
for i = 1:numGenesTreat
    fprintf(1,'%s, ',theGenesTreat{i});
end
fprintf(1,'\n');

%-------------------------------------------------------------------------------
% Now we look for differences in these columns compared to others:
isTreated = ismember(resultsTable.gene,theGenesTreat);
fprintf(1,'%u/%u ''treatment genes'' match GWAS-based analysis\n',...
                            sum(isTreated),numGenesTreat);

%-------------------------------------------------------------------------------
% Compare columns in resultsTable
propsToCompare = {'numGWASMapped','numLDSNPs','numPPIneighbors1DiseaseMapped',...
            'numPPIneighbors1DiseaseLD','percPPIneighbors1DiseaseMapped',...
            'percPPIneighbors1DiseaseLD','AllenMeanCoexp'};
numProps = length(propsToCompare);
f = figure('color','w');
for i = 1:numProps
    ax = subplot(2,4,i);
    dataVector = resultsTable.(propsToCompare{i});
    dataCell = cell(2,1);
    dataCell{1} = dataVector(isTreated);
    dataCell{2} = dataVector(~isTreated);
    [ff,xx] = BF_JitteredParallelScatter(dataCell,true,true,false);
    ax.XTick = 1:2;
    ax.XTickLabels = {sprintf('%streatment (%u)',whatDisease,sum(isTreated)),...
                    sprintf('not-treatment (%u)',sum(~isTreated))};
    title(propsToCompare{i});
end

end
