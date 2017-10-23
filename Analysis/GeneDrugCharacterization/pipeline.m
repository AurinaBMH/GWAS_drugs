% Pipeline for producing table characterizing individual genes

% Parameters:
PPINevidenceThreshold = 0; % evidence threshold for including PPI interactions
whatDisease = 'all'; % pick a disease to focus on: 'all','SZP','ASD','ADHD','BIP','MDD'

%===============================================================================
% LOAD AND PROCESS DATA from .csv files:
%===============================================================================

% List of all genes with drug targets in Drugbank are in **gene_ATC_matrix.csv**
[geneDrugTable,allUniqueGenes] = DrugGeneImport();
numUniqueGenes = length(allUniqueGenes);
numUniqueDrugs = length(unique(geneDrugTable.drugName));

% Import drug classification
drugClassTable = DrugClassImport();

% SNP, gene, disease, GWAS, LD annotations:
[SNPAnnotationTable,SNPGeneMap,allDiseaseSNPs] = SNPAnnotationImport(whatDisease);
fprintf(1,'%u/%u genes with drug targets have annotations\n',...
        sum(ismember(allUniqueGenes,SNPAnnotationTable.mappedGene)),numUniqueGenes);

isGWASAndNotEmpty = SNPAnnotationTable.isGWAS & cellfun(@(x)~isempty(x),SNPAnnotationTable.mappedGene);
allDiseaseGenesMapped = unique(SNPAnnotationTable.mappedGene(isGWASAndNotEmpty));
isLDAndNotEmpty = SNPAnnotationTable.isLD & cellfun(@(x)~isempty(x),SNPAnnotationTable.mappedGene);
allDiseaseGenesLD = unique(SNPAnnotationTable.mappedGene(isLDAndNotEmpty));

% Load the LD relationship data:
LDRelateTable = LDImport();

% Get cis-eQTL data:
[eQTLTable,isEGene,isSNPGene] = eQTLImport();

% Get PPIN data:
fileName = sprintf('PPI_Adj_th%u.mat',PPINevidenceThreshold);
try
    fprintf(1,'Loading PPIN data for evidence threshold of %u...',PPINevidenceThreshold);
    PPIN = load(fileName,'AdjPPI','geneNames');
    fprintf(1,' Loaded!\n');
catch
    error('No precomputed data for evidence threshold of %u... RECOMPUTING!!!',...
                    PPINevidenceThreshold)
    % PPIN = PPINImport(PPINevidenceThreshold);
end
% Load pairwise distances:
fileName = sprintf('PPI_Dist_th%u.mat',PPINevidenceThreshold);
try
    PPIND = load(fileName,'distMatrix');
    PPIN.distMatrix = PPIND.distMatrix;
    clear('PPIND');
catch
    warning('No PPI shortest path information found');
end
% Get indices of mapped disease genes on PPI network:
PPI_isDiseaseGeneMapped = cellfun(@(x)find(strcmp(x,PPIN.geneNames)),allDiseaseGenesMapped);

%-------------------------------------------------------------------------------
% Now loop through genes to characterize...
numGWASMapped = zeros(numUniqueGenes,1);
numLDSNPs = zeros(numUniqueGenes,1);
numEGenes = zeros(numUniqueGenes,1);
numSNPGenes = zeros(numUniqueGenes,1);
numEGenes_LD = zeros(numUniqueGenes,1);
numSNPGenes_LD = zeros(numUniqueGenes,1);
numLDeGeneseQTL = zeros(numUniqueGenes,1);
numLDSNPGeneseQTL = zeros(numUniqueGenes,1);
numPPIneighborsDiseaseMapped = zeros(numUniqueGenes,1);
numPPIneighborsDiseaseLD = zeros(numUniqueGenes,1);
meanPPIDistance = zeros(numUniqueGenes,1);
matchingDrugsString = cell(numUniqueGenes,1);

for i = 1:numUniqueGenes
    gene_i = allUniqueGenes{i};
    fprintf(1,'[%u/%u]: %s\n',i,numUniqueGenes,gene_i);

    %-------------------------------------------------------------------------------
    % Preliminaries:
    %-------------------------------------------------------------------------------
    % (I assume this can be comprehensive given the data provided??? Maybe not??)
    theLDgenes = GiveMeLDGenes(gene_i,SNPGeneMap,LDRelateTable,allDiseaseSNPs);
    % fprintf(1,'%u genes LD to the target\n',length(theLDgenes));

    %-------------------------------------------------------------------------------
    % --numGWASMapped: the number of GWAS SNPs mapped directly to the gene
    %-------------------------------------------------------------------------------
    isGeneAndMapped = strcmp(SNPAnnotationTable.mappedGene,gene_i) & SNPAnnotationTable.isGWAS;
    theMappedSNPs = unique(SNPAnnotationTable.SNP_id(isGeneAndMapped));
    if isempty(theMappedSNPs)
        numGWASMapped(i) = 0;
    else
        numGWASMapped(i) = length(theMappedSNPs);
    end

    %-------------------------------------------------------------------------------
    % --numLD: the number of GWAS SNPs LD to a gene
    %-------------------------------------------------------------------------------
    % If I want to count the number of SNPs that are LD to a given gene, I need to:
    % (i) Check every SNP in that gene (?? Using SNP_identifier.csv??)
    % (ii) Get a list of all (unique) SNPs are LD to these SNPs (Using LD_r2.csv)
    % (iii) Check all LD SNPs for direct disease annotations (in SNP_identifier.csv).
    isGeneAndLD = strcmp(SNPAnnotationTable.mappedGene,gene_i) & SNPAnnotationTable.isLD;
    theLDSNPs = unique(SNPAnnotationTable.SNP_id(isGeneAndLD));
    % NB: No way to double check that all LD annotations are from unique SNPs...
    % (is an assumption that Janette has formatted the file in this way)

    % Instead, we count all unique SNPs that are in a gene and LD to some GWAS hit:
    if isempty(theLDSNPs)
        numLDSNPs(i) = 0;
    else
        % Exclude SNPs that are mapped directly
        inBoth = intersect(theMappedSNPs,theLDSNPs);
        numLDSNPs(i) = length(theLDSNPs) - length(inBoth);
    end

    %-------------------------------------------------------------------------------
    % eQTLs-GWAS
    %-------------------------------------------------------------------------------
    isGene_i = strcmp(eQTLTable.eQTL,gene_i);
    % eGenes
    eGeneCandidates = unique(eQTLTable.SNP(isEGene & isGene_i));
    SNPGeneCandidates = unique(eQTLTable.SNP(isSNPGene & isGene_i));

    % Filter lists by disease (assume Janette already did??)
    isDisease_e = ismember(eGeneCandidates,allDiseaseSNPs);
    isDisease_SNP = ismember(SNPGeneCandidates,allDiseaseSNPs);

    if ~all(isDisease_e)
        warning('Some eGenes have no disease annotation??');
    end
    if ~all(isDisease_SNP)
        warning('Some SNPGenes have no disease annotation??');
    end

    % Count how many eGenes/SNPgenes have direct disease annotations
    numEGeneSNPs(i) = sum(isDisease_e);
    numSNPGeneSNPs(i) = sum(isDisease_SNP);
    % fprintf(1,'%u eGeneSNPS\n',numEGeneSNPs(i));
    % fprintf(1,'%u SNPGeneSNPS\n',numSNPGeneSNPs(i));

    %-------------------------------------------------------------------------------
    % eQTLs->SNP->LD
    %-------------------------------------------------------------------------------
    % How many unique SNPs are LD to SNPs that are eQTL for a given gene??
    eGeneSNPs = eGeneCandidates(isDisease_e);
    SNPGeneSNPs = SNPGeneCandidates(isDisease_SNP);

    % How many SNPs are LD to eGenes/SNPgenes
    % (1) Get all SNPs eQTL to the target gene
    % (2) count how many SNPs are LD to these eQTL SNPs
    theLDSNPs = cell(numEGeneSNPs(i),1);
    for j = 1:numEGeneSNPs(i)
        theLDSNPs{j} = GiveMeLDSNPs(eGeneSNPs{j},LDRelateTable,allDiseaseSNPs);
    end
    numEGenes_LD(i) = length(unique(vertcat(theLDSNPs{:})));

    % How many SNPs are LD to SNPgenes
    theLDSNPs = cell(numSNPGeneSNPs(i),1);
    for j = 1:numSNPGeneSNPs(i)
        theLDSNPs{j} = GiveMeLDSNPs(SNPGeneSNPs{j},LDRelateTable,allDiseaseSNPs);
    end
    numSNPGenes_LD(i) = length(unique(vertcat(theLDSNPs{:})));

    %-------------------------------------------------------------------------------
    % gene->LD->eQTL-SNPs
    %-------------------------------------------------------------------------------
    % How many SNPs are eQTL to SNPs LD to the gene of interest?
    % (1) get genes LD to the target gene
    % (2) count SNPs that are eQTL to LD genes
    numLDgenes = length(theLDgenes);
    eQTLSNPs_egene = cell(numLDgenes,1);
    eQTLSNPs_SNPgene = cell(numLDgenes,1);
    for j = 1:numLDgenes
        isGene_j = strcmp(eQTLTable.eQTL,theLDgenes{j});
        eGeneCandidates = unique(eQTLTable.SNP(isEGene & isGene_j));
        isDisease_e = ismember(eGeneCandidates,allDiseaseSNPs);
        SNPGeneCandidates = unique(eQTLTable.SNP(isSNPGene & isGene_j));
        isDisease_SNP = ismember(SNPGeneCandidates,allDiseaseSNPs);
        eQTLSNPs_egene{j} = eGeneCandidates(isDisease_e);
        eQTLSNPs_SNPgene{j} = SNPGeneCandidates(isDisease_SNP);
    end
    numLDeGeneseQTL(i) = length(unique(vertcat(eQTLSNPs_egene{:})));
    numLDSNPGeneseQTL(i) = length(unique(vertcat(eQTLSNPs_SNPgene{:})));

    %-------------------------------------------------------------------------------
    % PPIN Neighbors
    %-------------------------------------------------------------------------------
    % Genes that are neighbors on the PPI network:
    PPI_index = strcmp(PPIN.geneNames,gene_i);
    PPI_neighbors_index = PPIN.AdjPPI(PPI_index,:);
    PPI_neighbors_gene = PPIN.geneNames(PPI_neighbors_index);
    % How many first degree neighbors are in the disease list (mapped/LD)?:
    numPPIneighborsDiseaseMapped(i) = sum(ismember(PPI_neighbors_gene,allDiseaseGenesMapped));
    numPPIneighborsDiseaseLD(i) = sum(ismember(PPI_neighbors_gene,allDiseaseGenesLD));

    %-------------------------------------------------------------------------------
    % PPIN Distances
    %-------------------------------------------------------------------------------
    % What is the mean path length to genes on the disease list (mapped/LD)?:
    if isfield(PPIN,'distMatrix')
        allDistances = PPIN.distMatrix(PPI_index,PPI_isDiseaseGeneMapped);
        meanPPIDistance(i) = mean(allDistances);
    end

    %-------------------------------------------------------------------------------
    % Drugs, classes
    %-------------------------------------------------------------------------------
    % Match to drugs using geneDrugTable
    % assign class using drugClassTable
    matchingDrugs = geneDrugTable.drugName(strcmp(geneDrugTable.geneName,gene_i));
    assignClass = @(druggy)drugClassTable.whatClass(strcmp(drugClassTable.drugName,druggy));
    matchingDrugsClass = cellfun(@(x)sprintf('%s(%s)',x,char(assignClass(x))),...
                        matchingDrugs,'UniformOutput',false);
    matchingDrugsString{i} = BF_cat(matchingDrugsClass);
end

% Now make a table
gene = allUniqueGenes;
results = table(gene,numGWASMapped,numLDSNPs,numPPIneighborsDiseaseMapped,numPPIneighborsDiseaseLD,...
                meanPPIDistance,numEGenes,numSNPGenes,numEGenes_LD,numSNPGenes_LD,...
                numLDeGeneseQTL,numLDSNPGeneseQTL,matchingDrugsString);
% Sort by column, then by column, etc. in ordered hierarchy:
results = sortrows(results,{'numGWASMapped','numLDSNPs','numPPIneighborsDiseaseMapped',...
                'numPPIneighborsDiseaseLD','meanPPIDistance','numEGenes','numSNPGenes','numEGenes_LD',...
                'numSNPGenes_LD','numLDeGeneseQTL','numLDSNPGeneseQTL'},'descend');
fprintf(1,'%u genes have mapped and LD\n',sum(results.numLDSNPs > 0 & results.numGWASMapped > 0));
fprintf(1,'%u genes have mapped but no LD\n',sum(results.numLDSNPs==0 & results.numGWASMapped > 0));
fprintf(1,'%u genes have LD but no mapped\n',sum(results.numLDSNPs > 0 & results.numGWASMapped==0));
fprintf(1,'%u genes have eGene annotations\n',sum(results.numEGenes > 0));
fprintf(1,'%u genes have SNPgene annotations\n',sum(results.numSNPGenes > 0));
fprintf(1,'%u genes have eGene-LD annotations\n',sum(results.numEGenes_LD > 0));
fprintf(1,'%u genes have SNPgene-LD annotations\n',sum(results.numSNPGenes_LD > 0));
fprintf(1,'%u genes have LD-eGene annotations\n',sum(results.numEGenes_LD > 0));
fprintf(1,'%u genes have LD-SNPgene annotations\n',sum(results.numSNPGenes_LD > 0));

display(results(1:60,:));
