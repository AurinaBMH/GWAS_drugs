function [numPPIneighbors1,percPPIneighbors1,meanPPIDistance] = TellMePPIInfo(PPINevidenceThreshold,contextGenes,genesChar,doWeighted)
% Characterize each gene in the list genesChar in terms of its nearness to the
% set of context genes in contextGenes
%-------------------------------------------------------------------------------

if nargin < 4
    doWeighted = true;
end

%-------------------------------------------------------------------------------
% Filenames:
if doWeighted
    extraText = '_w';
    fprintf(1,'Loading weighted PPIN data...\n');
else
    extraText = sprintf('_th0%u',PPINevidenceThreshold*10);
    fprintf(1,'Loading PPIN data for evidence threshold of %.2f...\n',PPINevidenceThreshold);
end
fileNameAdj = sprintf('PPI_Adj%s.mat',extraText);
fileNameGeneLabels = sprintf('PPI_geneLabels%s.mat',extraText);
fileNamePDist = sprintf('PPI_Dist%s.mat',extraText);

%-------------------------------------------------------------------------------
% LOAD PPIN DATA (precomputed from PPINImport):
%-------------------------------------------------------------------------------
try
    load(fileNameAdj,'AdjPPI');
    load(fileNameGeneLabels,'geneNames');
    PPIN = struct();
    PPIN.AdjPPI = AdjPPI;
    PPIN.geneNames = geneNames; % (actually protein names)
    clear('AdjPPI','geneNames');
    fprintf(1,'Data loaded from %s and %s!\n',fileNameAdj,fileNameGeneLabels);
catch
    error('No precomputed data for evidence threshold of %.1f...!!!',...
                    PPINevidenceThreshold)
end

%-------------------------------------------------------------------------------
% 1. Map gene names (HGNC) -> names that match in the PPIN
% I've distinguished these as genesChar -> allUniqueProteins
%-------------------------------------------------------------------------------
% mapHow = 'UniProt'; % original method (no good)
mapHow = 'HGNC'; % Janette's new method for mapping... :)
switch mapHow
case 'UniProt'
    fprintf(1,'Mapping gene names to UniProt protein names. This actually makes things worse\n');
    [geneNameHGNC,proteinNameUniprot,allUniqueProteins] = ImportGeneUniProt(genesChar,PPIN.geneNames);
case 'HGNC'
    fprintf(1,'Mapping gene names to PPIN names (this only minimally helps)...\n');
    allUniqueProteins = ImportProteinGeneMapping(genesChar,PPIN.geneNames);
case 'exact'
    allUniqueProteins = genesChar;
end

%-------------------------------------------------------------------------------
% Get indices of mapped disease genes on PPI network:
%-------------------------------------------------------------------------------
PPI_isInContext = ismember(PPIN.geneNames,contextGenes);
fprintf(1,'%u/%u context genes are in the PPI network\n',sum(PPI_isInContext),length(contextGenes));
fprintf(1,'%u/%u genes to be characterized could be matched to PPI network data\n',...
                    sum(ismember(PPIN.geneNames,genesChar)),length(genesChar));

%-------------------------------------------------------------------------------
% Load shortest path distances on the PPI network:
% (cf. ComputePPIDist)
%-------------------------------------------------------------------------------
try
    PPIND = load(fileNamePDist,'distMatrix');
    PPIN.distMatrix = PPIND.distMatrix;
    clear('PPIND');
catch
    warning('No PPI shortest path information found');
end

%-------------------------------------------------------------------------------
% Initialize variables:
numGenesChar = length(genesChar);
numPPIneighbors1 = nan(numGenesChar,1);
percPPIneighbors1 = nan(numGenesChar,1);
meanPPIDistance = nan(numGenesChar,1);

%-------------------------------------------------------------------------------
% PPIN Neighbors
%-------------------------------------------------------------------------------
% Count disease genes that are 1-step neighbors on the PPI network:
% Match using "protein" nomenclature, mapped during processing steps above
for i = 1:numGenesChar
    gene_i = genesChar{i};
    protein_i = allUniqueProteins{i};
    PPI_index = find(strcmp(PPIN.geneNames,protein_i));

    %-------------------------------------------------------------------------------
    if isempty(PPI_index)
        fprintf('%s -> %s is not represented in the PPI network\n',gene_i,protein_i)
        % No match -- so cannot use any info from PPI network for this gene:
        continue;
    end

    %---------------------------------------------------------------------------
    % Get the 1-step neighbors
    PPI_neighbors_index = find(PPIN.AdjPPI(PPI_index,:));

    % Don't allow self-matches:
    PPI_neighbors_index(PPI_index) = 0;

    % Count neighbors:
    numNeighbors = length(PPI_neighbors_index);

    if numNeighbors==0
        numPPIneighbors1(i) = 0;
        percPPIneighbors1(i) = NaN;
    else
        % PPI_neighbors_index = union(PPI_neighbors_index,PPI_index); % include the target gene
        PPI_neighbors_gene = PPIN.geneNames(PPI_neighbors_index);
        % How many 1-step neighbors are in the disease list (mapped/LD)?:
        isInContext = ismember(PPI_neighbors_gene,contextGenes);
        numPPIneighbors1(i) = sum(isInContext);
        meanPPIneighbors1(i) = mean(isInContext);

        if doWeighted
            % What aggregrate evidence weight spread across 1-step neighbors are on the context list?:
            % (expressed as a proportion of total evidence weight):
            evidenceWeightContext = sum(PPIN.AdjPPI(PPI_index,isInContext));
            evidenceWeightTotal = sum(PPIN.AdjPPI(PPI_index,PPI_neighbors_index));
            percPPIneighbors1(i) = evidenceWeightContext/evidenceWeightTotal;
        else
            % As a percentage of the number of neighbors (~penalizes genes with many neighbors):
            percPPIneighbors1(i) = 100*meanPPIneighbors1;
            % fprintf(1,'%s: %u/%u neighbors from context\n',protein_i,numPPIneighbors1(i),numNeighbors);
        end
    end
    %-------------------------------------------------------------------------------
    % PPIN Distance:
    % What is the mean path length to genes on the disease list (mapped/LD)?:
    if isfield(PPIN,'distMatrix')
        allDistances = PPIN.distMatrix(PPI_index,PPI_isInContext);
        meanPPIDistance(i) = mean(allDistances);
    end
end

% %-------------------------------------------------------------------------------
% % Package up output nicely:
% outputVar = struct();
% outputVar.numPPIneighbors1 = numPPIneighbors1;
% outputVar.percPPIneighbors1 = percPPIneighbors1;
% outputVar.meanPPIDistance = meanPPIDistance;

%===============================================================================
% Genes that are 2-step neighbors on the PPI network -- doesn't make much sense (so many!):
% (e.g., in one example case there are 3235 1-step neighbors; 17992 2-step neighbors)
% (this basically gives coverage to the entire network in two steps...)
% neighborsIndex = find(PPI_neighbors_index); % list of indices
% neighborNeighbors = arrayfun(@(x)find(PPIN.AdjPPI(x,:)),neighborsIndex,'UniformOutput',false);
% twoStepNeighbors = unique(horzcat(neighborNeighbors{:}));
% PPI_neighbors_gene = PPIN.geneNames(PPI_neighbors_index);
% % How many 1/2-step neighbors are in the disease list (mapped/LD)?:
% numPPIneighbors2DiseaseMapped(i) = sum(ismember(PPI_neighbors_gene,contextGenes));
% numPPIneighbors2DiseaseLD(i) = sum(ismember(PPI_neighbors_gene,allDiseaseGenesLD));

end
