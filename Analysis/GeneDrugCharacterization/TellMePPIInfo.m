function geneStats = TellMePPIInfo(contextGenes,genesChar,doWeighted,PPINevidenceThreshold,numSteps)
% Characterize each gene in the list genesChar in terms of its nearness to the
% set of context genes in contextGenes
%-------------------------------------------------------------------------------
if nargin < 3
    doWeighted = true;
end
if nargin < 4
    PPINevidenceThreshold = 0;
end
if nargin < 5
    numSteps = 1;
end

whatPPIData = 'HGNC';
numContextGenes = length(contextGenes);

%-------------------------------------------------------------------------------
% Filenames:
if doWeighted
    extraText = '_w';
    fprintf(1,'Loading weighted PPIN data...\n');
else
    extraText = sprintf('_th%u',PPINevidenceThreshold);
    fprintf(1,'Loading PPIN data for evidence threshold of %.2f...\n',PPINevidenceThreshold);
end
switch whatPPIData
case 'HGNC'
    preText = 'PPI_HGNC';
end
fileNameAdj = sprintf('%s_Adj%s.mat',preText,extraText);
fileNameGeneLabels = sprintf('%s_geneLabels%s.mat',preText,extraText);
fileNamePDist = sprintf('%s_Dist%s.mat',preText,extraText);

%-------------------------------------------------------------------------------
% LOAD PPIN DATA (precomputed from PPINImport):
try
    load(fileNameAdj,'AdjPPI');
    load(fileNameGeneLabels,'geneNames');
    PPIN = struct();
    PPIN.AdjPPI = AdjPPI;
    PPIN.geneNames = geneNames; % (actually protein names)
    clear('AdjPPI','geneNames');
    fprintf(1,'Data loaded from %s and %s!\n',fileNameAdj,fileNameGeneLabels);
catch
    keyboard
    error('No precomputed data for specified PPI network...!!!')
end

%-------------------------------------------------------------------------------
% 1. Map gene names (HGNC) -> names that match in the PPIN
% I've distinguished these as genesChar -> allUniqueProteins
%-------------------------------------------------------------------------------
% mapHow = 'UniProt'; % original method (no good)
% mapHow = 'HGNC'; % Janette's new method for mapping... (needs proteins)
mapHow = 'exact';
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
fprintf(1,'%u/%u context genes are in the PPI network\n',sum(PPI_isInContext),numContextGenes);
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

numPPIneighbors = cell(numSteps,1);
percPPIneighbors = cell(numSteps,1);
hasPath = cell(numSteps,1);
meanPPIDistance = nan(numGenesChar,1);

% Initialize and precompute longer (binary) paths
for k = 1:numSteps
    numPPIneighbors{k} = nan(numGenesChar,1);
    percPPIneighbors{k} = nan(numGenesChar,1);
    if k==1
        hasPath{1} = PPIN.AdjPPI;
    else
        hasPath{k} = double(hasPath{1})^k; % matrix power
    end
    hasPath{k}(logical(eye(size(hasPath{k})))) = 0;
end

%-------------------------------------------------------------------------------
% PPIN Neighbors
%-------------------------------------------------------------------------------
% Count disease genes that are 1-step neighbors on the PPI network:
% Match using "protein" nomenclature, mapped during processing steps above
for i = 1:numGenesChar
    gene_i = genesChar{i};
    protein_i = allUniqueProteins{i};

    % First try matching the protein:
    PPI_index = find(strcmp(PPIN.geneNames,protein_i));
    % Then try matching the gene:
    if isempty(PPI_index)
        PPI_index = find(strcmp(PPIN.geneNames,gene_i));
    end

    %-------------------------------------------------------------------------------
    if isempty(PPI_index)
        fprintf('%s -> %s is not represented in the PPI network\n',gene_i,protein_i)
        continue;
    end

    %-------------------------------------------------------------------------------
    % Get the k-step neighbors:
    for k = 1:numSteps
        iskStepNeighbor = (hasPath{k}(PPI_index,:) > 0);
        numNeighbors = sum(iskStepNeighbor);
        if numNeighbors==0
            numPPIneighbors{k}(i) = 0;
            percPPIneighbors{k}(i) = NaN;
        else
            PPI_neighbors_gene = PPIN.geneNames(iskStepNeighbor);
            isInContext = ismember(PPI_neighbors_gene,contextGenes);

            % How many k-step neighbors are in the context/disease list?:
            numPPIneighbors{k}(i) = sum(isInContext);

            if doWeighted
                % What aggregrate evidence weight spread across 1-step neighbors are on the context list?:
                % (expressed as a proportion of total evidence weight):
                evidenceWeightContext = sum(hasPath{k}(PPI_index,isInContext));
                evidenceWeightTotal = sum(hasPath{k}(PPI_index,isNeighbor));
                percPPIneighbors{k}(i) = evidenceWeightContext/evidenceWeightTotal;
            else
                % As a percentage of the number of neighbors (~penalizes genes with many neighbors):
                percPPIneighbors{k}(i) = 100*mean(isInContext);
                % fprintf(1,'%s: %u/%u neighbors from context\n',protein_i,numPPIneighbors1(i),numNeighbors);
            end
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

%===============================================================================
% Genes that are 2-step neighbors on the PPI network -- doesn't make much sense (so many!):
% (e.g., in one example case there are 3235 1-step neighbors; 17992 2-step neighbors)
% (this basically gives coverage to the entire network in two steps...)

%-------------------------------------------------------------------------------
% Package up output nicely:
geneStats = struct();
for k = 1:numSteps
    geneStats.sprintf('numPPIneighbors%u',k) = numPPIneighbors{k};
    geneStats.sprintf('percPPIneighbors%u',k) = percPPIneighbors{k};
end
geneStats.meanPPIDistance = meanPPIDistance;
if numSteps > 1
    geneStats.numPPIneighbors2 = numPPIneighbors2;
    geneStats.percPPIneighbors2 = percPPIneighbors2;
end

end
