function geneStats = TellMePPIInfo(mappedGenes,targetGenes,doWeighted,evidenceThreshold,numSteps)
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
numContextGenes = length(mappedGenes);

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

switch whatPPIData
case 'HGNCmatch'
    % No need to map out to other nomenclatures because data is already in the format
    allUniqueProteins = targetGenes;
otherwise
    warning('Unknown dataset, how should I match gene names -> protein names for optimal matching?')
end


%-------------------------------------------------------------------------------
% Get indices of mapped disease genes on PPI network:
%-------------------------------------------------------------------------------
PPI_isInContext = ismember(mappedGenes,PPIN.geneNames);
PPI_isInChar = ismember(targetGenes,PPIN.geneNames);
missingContext = find(~PPI_isInContext);
missingChar = find(~PPI_isInChar);
fprintf(1,'%u/%u context genes are in the PPI network\n',sum(PPI_isInContext),numContextGenes);
fprintf(1,'%u/%u genes to be characterized could be matched to PPI network data\n',...
                    sum(PPI_isInChar),length(targetGenes));

for i = 1:length(missingContext)
    fprintf(1,'[%u/%u context]: %s missing in PPIN\n',i,length(missingContext),...
                                mappedGenes{missingContext(i)});
end
for i = 1:length(missingChar)
    fprintf(1,'[%u/%u characterize]: %s missing in PPIN\n',i,length(missingChar),...
                                targetGenes{missingChar(i)});
end

%-------------------------------------------------------------------------------
% Initialize variables:
numtargetGenes = length(targetGenes);
numPPIneighbors = cell(numSteps,1);
gwasPPIneighbors = cell(numSteps,1);
percPPIneighbors = cell(numSteps,1);
weiPPIneighbors = cell(numSteps,1);
numCOMMONneighbors = cell(numSteps,1);
percCOMMONneighbors = cell(numSteps,1);
expWeiPPIneighbors = cell(numSteps,1);
weigwasPPIneighbors = cell(numSteps,1);
expWeigwasPPIneighbors = cell(numSteps,1);
hasPath = cell(numSteps,1);
meanPPIDistance = nan(numtargetGenes,1);
medianPPIDistance = nan(numtargetGenes,1);

% Initialize longer (binary) paths
for k = 1:numSteps
    numPPIneighbors{k} = nan(numtargetGenes,1);
    percPPIneighbors{k} = nan(numtargetGenes,1);
    weiPPIneighbors{k} = zeros(numtargetGenes,1);
    expWeiPPIneighbors{k} = zeros(numtargetGenes,1);
    numCOMMONneighbors{k} = nan(numtargetGenes,1);
    percCOMMONneighbors{k} = nan(numtargetGenes,1);
    weigwasPPIneighbors{k} = zeros(numtargetGenes,1);
    expWeigwasPPIneighbors{k} = zeros(numtargetGenes,1);
    gwasPPIneighbors{k} = nan(numtargetGenes,1);
end

%-------------------------------------------------------------------------------
% PPIN Neighbors
%-------------------------------------------------------------------------------
% Count disease genes that are 1-step neighbors on the PPI network:
% Match using "protein" nomenclature, mapped during processing steps above
fprintf(1,'Characterizing the properties of %u difference genes in the context of %u genes\n',...
                        numtargetGenes,numContextGenes);
% for a selected GWAS list, find a list of k-step neighbors
PPI_neighborsK = find_PPIneighbors(PPIN, mappedGenes, numSteps); 

for i = 1:numtargetGenes
    gene_i = targetGenes{i};
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
        
        % Sum of the number of labeled genes that are within k steps (maximum score is the number of labeled genes)

        numNeighbors = sum(iskStepNeighbor);
        if numNeighbors==0
            numPPIneighbors{k}(i) = 0;
            percPPIneighbors{k}(i) = NaN;
        else
            PPI_neighbors_gene = PPIN.geneNames(iskStepNeighbor);
            isInContext = ismember(PPI_neighbors_gene,mappedGenes);
            % find common neighbors instead of GWAS hits, use a list of the
            % neighbors of GWAS hits at a selected k threshold PPI_neighborsK{k}
            isInContextNeighbors = ismember(PPI_neighbors_gene,PPI_neighborsK{k});
            
            % How many k-step neighbors are in the context/disease list?:
            numPPIneighbors{k}(i) = sum(isInContext);
            numCOMMONneighbors{k}(i) = sum(isInContextNeighbors);

            if doWeighted
                % What aggregrate evidence weight spread across 1-step neighbors are on the context list?:
                % (expressed as a proportion of total evidence weight):
                evidenceWeightContext = sum(hasPath{k}(PPI_index,isInContext));
                evidenceWeightTotal = sum(hasPath{k}(PPI_index,isNeighbor));
                percPPIneighbors{k}(i) = evidenceWeightContext/evidenceWeightTotal;

            else
                % As a percentage of the number of neighbors (~penalizes genes with many neighbors):
                percPPIneighbors{k}(i) = 100*mean(isInContext);
                % for common neighbors use the union of both lists
                all_neighbors = length(unique(vertcat(PPI_neighbors_gene,PPI_neighborsK{k})));
                percCOMMONneighbors{k}(i) = 100*(sum(isInContextNeighbors)/all_neighbors); 
                % how many labeled (fromGWAS) neighbors a gene has / # of all labeled genes
                gwasPPIneighbors{k}(i) = sum(isInContext)/numContextGenes;
                % fprintf(1,'%s: %u/%u neighbors from context\n',protein_i,numPPIneighbors1(i),numNeighbors);
                
            end
        end

        % (ii) Weight contribution of neighbors by their pathlength:
        for kk = 1:k
            iskkStepNeighbor = (PPIN.distMatrix(PPI_index,:) == k);
            prop_kk_neighbors = 100*mean(ismember(PPIN.geneNames(iskkStepNeighbor),mappedGenes));
            propgwas_kk_neighbors = sum(ismember(PPIN.geneNames(iskkStepNeighbor),mappedGenes))/numContextGenes; 
            if ~isnan(prop_kk_neighbors)
                weiPPIneighbors{k}(i) = weiPPIneighbors{k}(i) + prop_kk_neighbors/kk;
                expWeiPPIneighbors{k}(i) = expWeiPPIneighbors{k}(i) + prop_kk_neighbors/factorial(kk);
                weigwasPPIneighbors{k}(i) = weigwasPPIneighbors{k}(i) + propgwas_kk_neighbors/kk; 
                expWeigwasPPIneighbors{k}(i) = weigwasPPIneighbors{k}(i) + propgwas_kk_neighbors/factorial(kk);
                                
                %Weighted sum of the distance to each labeled gene: e.g., 
                %if gene 1 is 1 step away, it gets scored 1 from this, 
                %then gene 2 is 3 steps away, it gets scored 1/3! for this gene, then gene 3, 4, ? 
                %So the maximum score is all genes are 1 step away. Can normalize by this.
            end
        end
    end

    %-------------------------------------------------------------------------------
    % PPIN Distance:
    % What is the mean path length to genes on the disease list (mapped/LD)?:
    isInContext = ismember(PPIN.geneNames,mappedGenes);
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
    
    geneStats.(sprintf('numCOMMONneighbors%u',k)) = numCOMMONneighbors{k};
    geneStats.(sprintf('percCOMMONneighbors%u',k)) = percCOMMONneighbors{k};
    
    geneStats.(sprintf('gwasPPIneighbors%u',k)) = gwasPPIneighbors{k};
    geneStats.(sprintf('weigwasPPIneighbors%u',k)) = weigwasPPIneighbors{k};
    geneStats.(sprintf('expWeigwasPPIneighbors%u',k)) = expWeigwasPPIneighbors{k};
    geneStats.(sprintf('gwasPPIneighbors%u',k)) = gwasPPIneighbors{k};

end
geneStats.medianPPIDistance = medianPPIDistance;
geneStats.meanPPIDistance = meanPPIDistance;

end
