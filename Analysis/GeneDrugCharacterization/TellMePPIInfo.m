function geneStats = TellMePPIInfo(contextGenes,genesChar,doWeighted,evidenceThreshold,numSteps)
% Characterize each gene in the list genesChar in terms of its nearness to the
% set of context genes in contextGenes
%-------------------------------------------------------------------------------
if nargin < 3
    doWeighted = false;
end
if nargin < 4
    evidenceThreshold = 0;
end
if nargin < 5
    numSteps = 5;
    fprintf(1,'%u steps by default\n',numSteps);
end

whatPPIData = 'HGNCmatch';
numContextGenes = length(contextGenes);

%-------------------------------------------------------------------------------
% Filenames:
fileNames = PPIFileNames(doWeighted,evidenceThreshold,whatPPIData);
fileNameGeneLabels = fileNames{2};
fileNameAdj = fileNames{3};
fileNamePDist = fileNames{4};

%-------------------------------------------------------------------------------
% LOAD PPIN DATA (precomputed from PPINImport):
%-------------------------------------------------------------------------------
try
    load(fileNameAdj,'AdjPPI');
    load(fileNameGeneLabels,'geneNames');
    PPIN = struct();
    PPIN.AdjPPI = AdjPPI;
    PPIN.geneNames = geneNames; % (/protein names?)
    clear('AdjPPI','geneNames');
    fprintf(1,'Data loaded from %s and %s!\n',fileNameAdj,fileNameGeneLabels);
catch
    error('No precomputed data for specified PPI network...!!!')
end
% Load shortest path distances on the PPI network (cf. ComputePPIDist):
try
    load(fileNamePDist,'distMatrix');
    PPIN.distMatrix = distMatrix;
    clear('distMatrix');
catch
    error('No PPI shortest path information found in %s, you must run ComputePPIDist',fileNamePDist);
end

%-------------------------------------------------------------------------------
% 1. Map gene names (HGNC) -> names that match in the PPIN
% I've distinguished these as genesChar -> allUniqueProteins
%-------------------------------------------------------------------------------
% mapHow = 'UniProt'; % original method (no good)
% mapHow = 'HGNC_map'; % Janette's new method for mapping... (needs proteins)
% mapHow = 'exact';
switch whatPPIData
case 'HGNCmatch'
    % No need to map out to other nomenclatures because data is already in the format
    allUniqueProteins = genesChar;
otherwise
    warning('Unknown dataset, how should I match gene names -> protein names for optimal matching?')
end
% switch mapHow
% case 'UniProt'
%     fprintf(1,'Mapping gene names to UniProt protein names. This actually makes things worse\n');
%     [geneNameHGNC,proteinNameUniprot,allUniqueProteins] = ImportGeneUniProt(genesChar,PPIN.geneNames);
% case 'HGNC_map'
%     fprintf(1,'Mapping gene names to PPIN names (this only minimally helps)...\n');
%     allUniqueProteins = ImportProteinGeneMapping(genesChar,PPIN.geneNames);
% end

%-------------------------------------------------------------------------------
% Get indices of mapped disease genes on PPI network:
%-------------------------------------------------------------------------------
PPI_isInContext = ismember(contextGenes,PPIN.geneNames);
PPI_isInChar = ismember(genesChar,PPIN.geneNames);
missingContext = find(~PPI_isInContext);
missingChar = find(~PPI_isInChar);
fprintf(1,'%u/%u context genes are in the PPI network\n',sum(PPI_isInContext),numContextGenes);
fprintf(1,'%u/%u genes to be characterized could be matched to PPI network data\n',...
                    sum(PPI_isInChar),length(genesChar));

for i = 1:length(missingContext)
    fprintf(1,'[%u/%u context]: %s missing in PPIN\n',i,length(missingContext),...
                                contextGenes{missingContext(i)});
end
for i = 1:length(missingChar)
    fprintf(1,'[%u/%u characterize]: %s missing in PPIN\n',i,length(missingChar),...
                                genesChar{missingChar(i)});
end

%-------------------------------------------------------------------------------
% Initialize variables:
numGenesChar = length(genesChar);
numPPIneighbors = cell(numSteps,1);
percPPIneighbors = cell(numSteps,1);
hasPath = cell(numSteps,1);
meanPPIDistance = nan(numGenesChar,1);
medianPPIDistance = nan(numGenesChar,1);

% Initialize longer (binary) paths
for k = 1:numSteps
    numPPIneighbors{k} = nan(numGenesChar,1);
    percPPIneighbors{k} = nan(numGenesChar,1);
    weiPPIneighbors{k} = zeros(numGenesChar,1);
    expWeiPPIneighbors{k} = zeros(numGenesChar,1);
end

%-------------------------------------------------------------------------------
% PPIN Neighbors
%-------------------------------------------------------------------------------
% Count disease genes that are 1-step neighbors on the PPI network:
% Match using "protein" nomenclature, mapped during processing steps above
fprintf(1,'Characterizing the properties of %u difference genes in the context of %u genes\n',...
                        numGenesChar,numContextGenes);
for i = 1:numGenesChar
    gene_i = genesChar{i};
    protein_i = allUniqueProteins{i};

    % First try matching the gene:
    PPI_index = find(strcmp(PPIN.geneNames,gene_i));
    % Then try matching the protein:
    if isempty(PPI_index)
        PPI_index = find(strcmp(PPIN.geneNames,protein_i));
    end
    if isempty(PPI_index)
        fprintf('%s/%s is not represented in the PPI network\n',gene_i,protein_i)
        continue;
    end

    %-------------------------------------------------------------------------------
    % Get the k-step neighbors:
    for k = 1:numSteps
        % (i) Don't distinguish neighbors by pathlength (take union):
        iskStepNeighbor = (PPIN.distMatrix(PPI_index,:) <= k);

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

        % (ii) Weight contribution of neighbors by their pathlength:
        for kk = 1:k
            iskkStepNeighbor = (PPIN.distMatrix(PPI_index,:) == k);
            prop_kk_neighbors = 100*mean(ismember(PPIN.geneNames(iskkStepNeighbor),contextGenes));
            if ~isnan(prop_kk_neighbors)
                weiPPIneighbors{k}(i) = weiPPIneighbors{k}(i) + prop_kk_neighbors/kk;
                expWeiPPIneighbors{k}(i) = expWeiPPIneighbors{k}(i) + prop_kk_neighbors/factorial(kk);
            end
        end
    end

    %-------------------------------------------------------------------------------
    % PPIN Distance:
    % What is the mean path length to genes on the disease list (mapped/LD)?:
    isInContext = ismember(PPIN.geneNames,contextGenes);
    allDistances = PPIN.distMatrix(PPI_index,isInContext);
    meanPPIDistance(i) = nanmean(allDistances);
    medianPPIDistance(i) = nanmedian(allDistances);
end

%-------------------------------------------------------------------------------
% Package up output nicely:
geneStats = struct();
for k = 1:numSteps
    geneStats.(sprintf('numPPIneighbors%u',k)) = numPPIneighbors{k};
    geneStats.(sprintf('percPPIneighbors%u',k)) = percPPIneighbors{k};
    geneStats.(sprintf('weiPPIneighbors%u',k)) = weiPPIneighbors{k};
    geneStats.(sprintf('expWeiPPIneighbors%u',k)) = expWeiPPIneighbors{k};
end
geneStats.medianPPIDistance = medianPPIDistance;
geneStats.meanPPIDistance = meanPPIDistance;

end
