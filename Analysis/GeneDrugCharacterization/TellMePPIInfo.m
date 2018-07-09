function [numPPIneighbors1,percPPIneighbors1,meanPPIDistance] = TellMePPIInfo(PPINevidenceThreshold,contextGenes,genesChar)
% Characterize each gene in the list genesChar in terms of its nearness to the
% set of context genes in contextGenes
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% LOAD PPIN DATA (precomputed from PPINImport):
%-------------------------------------------------------------------------------
fileNameAdj = sprintf('PPI_Adj_th0%u.mat',PPINevidenceThreshold*10);
fileNameGeneLabels = sprintf('PPI_geneLabels_th0%u.mat',PPINevidenceThreshold*10);
try
    fprintf(1,'Loading PPIN data for evidence threshold of %.2f...',PPINevidenceThreshold);
    load(fileNameAdj,'AdjPPI');
    load(fileNameGeneLabels,'geneNames');
    PPIN = struct();
    PPIN.AdjPPI = AdjPPI;
    PPIN.geneNames = geneNames; % (actually protein names)
    clear('AdjPPI','geneNames');
    fprintf(1,' Loaded from %s and %s!\n',fileNameAdj,fileNameGeneLabels);
catch
    error('No precomputed data for evidence threshold of %.1f...!!!',...
                    PPINevidenceThreshold)
    % PPIN = PPINImport(PPINevidenceThreshold);
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
    fprintf(1,'Mapping gene names to PPIN names. This only minimally helps...\n');
    allUniqueProteins = ImportProteinGeneMapping(genesChar,PPIN.geneNames);
case 'exact'
    allUniqueProteins = genesChar;
end

%-------------------------------------------------------------------------------
% Get indices of mapped disease genes on PPI network:
PPI_isDiseaseGeneMapped = ismember(PPIN.geneNames,contextGenes);
fprintf(1,'%u/%u GWAS disease genes are in the PPI network\n',...
            sum(PPI_isDiseaseGeneMapped),length(contextGenes));
PPI_isDrugGene = ismember(PPIN.geneNames,genesChar);
fprintf(1,'%u/%u drug-target genes could be matched to PPI network data\n',sum(PPI_isDrugGene),...
                                        length(genesChar));

%-------------------------------------------------------------------------------
% Load shortest path distances on the PPI network:
fileNamePDist = sprintf('PPI_Dist_th0%u.mat',PPINevidenceThreshold*10);
try
    PPIND = load(fileNamePDist,'distMatrix');
    PPIN.distMatrix = PPIND.distMatrix;
    clear('PPIND');
catch
    warning('No PPI shortest path information found');
end
% list the bad ones:
% noGood = find(~isPPIGene);
% for i = 1:length(noGood)
%     fprintf(1,'%s\n',genesChar{noGood(i)});
% end

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
        warning('%s->%s not represented in the PPI network',gene_i,protein_i)
        % No match -- so cannot use any info from PPI network for this gene:
        continue;
    end

    % Get the 1-step neighbors
    PPI_neighbors_index = find(PPIN.AdjPPI(PPI_index,:));
    numNeighbors = length(PPI_neighbors_index);
    if numNeighbors==0
        numPPIneighbors1DiseaseMapped(i) = 0;
        numPPIneighbors1DiseaseLD(i) = 0;
        percPPIneighbors1DiseaseMapped(i) = NaN;
        percPPIneighbors1DiseaseLD(i) = NaN;
    else
        % PPI_neighbors_index = union(PPI_neighbors_index,PPI_index); % include the target gene
        PPI_neighbors_gene = PPIN.geneNames(PPI_neighbors_index);

        % How many 1-step neighbors are in the disease list (mapped/LD)?:
        numPPIneighbors1(i) = sum(ismember(PPI_neighbors_gene,contextGenes));

        % As a percentage of the number of neighbors (~penalizes genes with many neighbors):
        percPPIneighbors1(i) = 100*numPPIneighbors1(i)/length(PPI_neighbors_index);

        fprintf(1,'%s: %u/%u neighbors from context\n',protein_i,numPPIneighbors1(i),numNeighbors);
    end
    %-------------------------------------------------------------------------------
    % PPIN Distance:
    % What is the mean path length to genes on the disease list (mapped/LD)?:
    if isfield(PPIN,'distMatrix')
        allDistances = PPIN.distMatrix(PPI_index,PPI_isDiseaseGeneMapped);
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
