% Pipeline for producing table characterizing individual genes

%===============================================================================
% LOAD AND PROCESS DATA from .csv files:
%===============================================================================

% List of all genes with drug targets in Drugbank are in **gene_ATC_matrix.csv**
[geneDrugTable,allUniqueGenes] = DrugGeneImport();
numUniqueGenes = length(allUniqueGenes);
numUniqueDrugs = length(unique(geneDrug.drugName));

% SNP, gene, disease, GWAS, LD annotations:
[SNPAnnotationTable,SNPGeneMap,allDiseaseSNPs] = SNPAnnotationImport()

% Load the LD relationship data:
LDRelateTable = LDImport();

% Get cis-eQTL data:
[eQTLTable,isEGene,isSNPGene] = eQTLImport();

%-------------------------------------------------------------------------------
% Infer the LD gene (i.e., the gene causing the annotation) for LD annotations
% in SNPAnnotationTable
%-------------------------------------------------------------------------------
% For LD SNPs, we can annotate the SNP that caused the annotation
% LDGene = cell(numAnnotations,1);
% for i = 1:numAnnotations
%     if SNPAnnotationTable.isLD(i)
%         % We want to determine which gene caused the annotation:
%         theLDSNP1 = LDRelateTable.SNP_id_2(strcmp(LDRelateTable.SNP_id_1,SNPAnnotationTable.SNP_id{i}));
%         theLDSNP2 = LDRelateTable.SNP_id_1(strcmp(LDRelateTable.SNP_id_2,SNPAnnotationTable.SNP_id{i}));
%         theLDSNPs = unique(vertcat(theLDSNP1,theLDSNP2));
%         if ischar(theLDSNPs), theLDSNPs = {theLDSNPs}; end
%
%         % Match to genes:
%         theLDGenes = cellfun(@(x)SNPGeneMap.mappedGene(strcmp(SNPGeneMap.SNP_id,x)),theLDSNPs,'UniformOutput',false);
%
%         if isempty(theLDGenes)
%             % No matches? Must be an LD relationship to a SNP without a gene
%             warning('No LD gene for %s',SNPAnnotationTable.SNP_id{i});
%             LDGene{i} = '';
%         elseif length(theLDGenes)==1
%             % One match -- easy:
%             LDGene{i} = theLDGenes{1};
%         else
%             fprintf(1,'Multiple matches :-O\n');
%         end
%
%     else
%         LDGene{i} = '';
%     end
% end

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

for i = 1:numUniqueGenes
    gene_i = allUniqueGenes{i};
    fprintf(1,'[%u/%u]: %s\n',i,numUniqueGenes,gene_i);

    %-------------------------------------------------------------------------------
    % Preliminaries:
    %-------------------------------------------------------------------------------
    % (I assume this can be comprehensive given the data provided??? Maybe not??)
    theLDgenes = GiveMeLDGenes(gene_i,SNPGeneMap,LDRelateTable,allDiseaseSNPs);
    fprintf(1,'%u genes LD to the target\n',length(theLDgenes));

    %-------------------------------------------------------------------------------
    % --numGWASMapped: the number of GWAS SNPs mapped directly to a gene
    %-------------------------------------------------------------------------------
    theSNPs = SNPAnnotationTable.SNP_id(strcmp(SNPAnnotationTable.mappedGene,gene_i) & SNPAnnotationTable.isGWAS);
    if isempty(theSNPs)
        numGWASMapped(i) = 0;
    else
        numGWASMapped(i) = length(unique(theSNPs));
    end

    %-------------------------------------------------------------------------------
    % --numLD: the number of GWAS SNPs LD to a gene
    %-------------------------------------------------------------------------------
    % If I want to count the number of SNPs that are LD to a given gene, I need to:
    % (i) Check every SNP in that gene (?? Using SNP_identifier.csv??)
    % (ii) Get a list of all (unique) SNPs are LD to these SNPs (Using LD_r2.csv)
    % (iii) Check all LD SNPs for direct disease annotations (in SNP_identifier.csv).
    isGeneAndLD = strcmp(SNPAnnotationTable.mappedGene,gene_i) & SNPAnnotationTable.isLD;
    % NB: No way to double check that all LD annotations are from unique SNPs...
    % (is an assumption that Janette has formatted the file in this way)
    numLDSNPs(i) = sum(isGeneAndLD);

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
    % PPIN
    %-------------------------------------------------------------------------------

end

% Now make a table
gene = allUniqueGenes;
results = table(gene,numGWASMapped,numLDSNPs,numEGenes,numSNPGenes,...
                numEGenes_LD,numSNPGenes_LD,numLDeGeneseQTL,numLDSNPGeneseQTL);
results = sortrows(results,{'numGWASMapped','numLDSNPs','numEGenes','numSNPGenes',...
        'numEGenes_LD','numSNPGenes_LD','numLDeGeneseQTL','numLDSNPGeneseQTL'},'descend');
fprintf(1,'%u genes have mapped and LD\n',sum(results.numLDSNPs > 0 & results.numGWASMapped > 0));
fprintf(1,'%u genes have mapped but no LD\n',sum(results.numLDSNPs==0 & results.numGWASMapped > 0));
fprintf(1,'%u genes have LD but no mapped\n',sum(results.numLDSNPs > 0 & results.numGWASMapped==0));
fprintf(1,'%u genes have eGene annotations\n',sum(results.numEGenes > 0));
fprintf(1,'%u genes have SNPgene annotations\n',sum(results.numSNPGenes > 0));
fprintf(1,'%u genes have eGene-LD annotations\n',sum(results.numEGenes_LD > 0));
fprintf(1,'%u genes have SNPgene-LD annotations\n',sum(results.numSNPGenes_LD > 0));
fprintf(1,'%u genes have LD-eGene annotations\n',sum(results.numEGenes_LD > 0));
fprintf(1,'%u genes have LD-SNPgene annotations\n',sum(results.numSNPGenes_LD > 0));
