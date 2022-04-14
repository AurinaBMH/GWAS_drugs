function outputVar = TellMeQTL(allDiseaseGenesMapped,allDiseaseSNPs,allUniqueGenes)

%-------------------------------------------------------------------------------
% Get cis-eQTL data:
[eQTLTable,isEGene,isSNPGene] = eQTLImport();

%-------------------------------------------------------------------------------
% Load the LD relationship data:
LDRelateTable = LDImport();

%-------------------------------------------------------------------------------
% Initialize variables:
numEGenes = nan(numUniqueGenes,1);
numSNPGenes = nan(numUniqueGenes,1);
numEGenes_LD = nan(numUniqueGenes,1);
numSNPGenes_LD = nan(numUniqueGenes,1);
numLDeGeneseQTL = nan(numUniqueGenes,1);
numLDSNPGeneseQTL = nan(numUniqueGenes,1);

%-------------------------------------------------------------------------------
% Start looking across allUniqueGenes
%-------------------------------------------------------------------------------
numUniqueGenes = length(allUniqueGenes);
for i = 1:numUniqueGenes
    gene_i = allUniqueGenes{i};

    %-------------------------------------------------------------------------------
    % eQTL -> GWAS
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
    % eQTLs -> SNP -> LD
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
    % gene -> LD -> eQTL -> SNPs
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

end

%-------------------------------------------------------------------------------
% Package output up nicely:
outputVar = struct();
outputVar.numEGenes = numEGenes;
outputVar.numSNPGenes = numSNPGenes;
outputVar.numEGenes_LD = numEGenes_LD;
outputVar.numSNPGenes_LD = numSNPGenes_LD;
outputVar.numLDeGeneseQTL = numLDeGeneseQTL;
outputVar.numLDSNPGeneseQTL = numLDSNPGeneseQTL;

end
