function [numGWASMapped,numLDSNPs] = TellMeDNADistance(genesChar,SNPAnnotationTable)
% Looks at direct mappings and LD-vicinity of a set of genes of interest (genesChar)
% using gene information from SNPAnnotationTable
%-------------------------------------------------------------------------------

numGenesChar = length(genesChar);
%-------------------------------------------------------------------------------
% Initialize:
numGWASMapped = zeros(numGenesChar,1);
numLDSNPs = zeros(numGenesChar,1);

for i = 1:numGenesChar
    gene_i = genesChar{i};

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
end

end
