function geneStats = TellMePPIInfo(contextGenes,genesChar,doWeighted,PPINevidenceThreshold)
% Characterize each gene in the list genesChar in terms of its nearness to the
% set of context genes in contextGenes
%-------------------------------------------------------------------------------
if nargin < 3
    doWeighted = true;
end
if nargin < 4
    PPINevidenceThreshold = 0;
end

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
    error('No precomputed data for specified PPI network...!!!')
end

%-------------------------------------------------------------------------------
% 1. Map gene names (HGNC) -> names that match in the PPIN
% I've distinguished these as genesChar -> allUniqueProteins
%-------------------------------------------------------------------------------
% mapHow = 'UniProt'; % original method (no good)
mapHow = 'HGNC'; % Janette's new method for mapping... (needs proteins)
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
numPPIneighbors1 = nan(numGenesChar,1);
percPPIneighbors1 = nan(numGenesChar,1);
meanPPIDistance = nan(numGenesChar,1);
if ~doWeighted
    numPPIneighbors2 = nan(numGenesChar,1);
    percPPIneighbors2 = nan(numGenesChar,1);
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

    %---------------------------------------------------------------------------
    % Get the 1-step neighbors
    isNeighbor = (PPIN.AdjPPI(PPI_index,:) > 0);
    % Don't allow self-matches:
    isNeighbor(PPI_index) = 0;
    neighborsIndex = find(isNeighbor);

    % Count neighbors:
    numNeighbors = length(neighborsIndex);

    if numNeighbors==0
        numPPIneighbors1(i) = 0;
        percPPIneighbors1(i) = NaN;
    else
        % neighborsIndex = union(neighborsIndex,PPI_index); % include the target gene
        PPI_neighbors_gene = PPIN.geneNames(neighborsIndex);
        % How many 1-step neighbors are in the disease list (mapped/LD)?:
        isInContext = ismember(PPI_neighbors_gene,contextGenes);
        numPPIneighbors1(i) = sum(isInContext);

        if doWeighted
            % What aggregrate evidence weight spread across 1-step neighbors are on the context list?:
            % (expressed as a proportion of total evidence weight):
            evidenceWeightContext = sum(PPIN.AdjPPI(PPI_index,isInContext));
            evidenceWeightTotal = sum(PPIN.AdjPPI(PPI_index,neighborsIndex));
            percPPIneighbors1(i) = evidenceWeightContext/evidenceWeightTotal;
        else
            % As a percentage of the number of neighbors (~penalizes genes with many neighbors):
            percPPIneighbors1(i) = 100*mean(isInContext);
            % fprintf(1,'%s: %u/%u neighbors from context\n',protein_i,numPPIneighbors1(i),numNeighbors);
        end
    end

    %-------------------------------------------------------------------------------
    % Get the (up to) two-step neighbors
    if ~doWeighted
        neighborNeighbors = arrayfun(@(x)find(PPIN.AdjPPI(x,:)),neighborsIndex,'UniformOutput',false);
        twoStepNeighbors = unique(horzcat(neighborNeighbors{:}));
        numNeighbors = length(twoStepNeighbors);
        if numNeighbors==0
            numPPIneighbors2(i) = 0;
            percPPIneighbors2(i) = NaN;
        else
            gene_2step_neighbors = PPIN.geneNames(twoStepNeighbors);
            isInContext = ismember(gene_2step_neighbors,contextGenes);
            % How many 2-step neighbors are in the context:
            numPPIneighbors2(i) = sum(isInContext);
            percPPIneighbors2(i) = 100*mean(isInContext);
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
geneStats.numPPIneighbors1 = numPPIneighbors1;
geneStats.percPPIneighbors1 = percPPIneighbors1;
geneStats.meanPPIDistance = meanPPIDistance;
if ~doWeighted
    geneStats.numPPIneighbors2 = numPPIneighbors2;
    geneStats.percPPIneighbors2 = percPPIneighbors2;
end

end
