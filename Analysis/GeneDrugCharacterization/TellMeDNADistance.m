function [numGWASMapped,numLDSNPs] = TellMeDNADistance(genesChar,SNPAnnotationTable,LDthreshold)
% Looks at direct mappings and LD-vicinity of a set of genes of interest (genesChar)
% using gene information from SNPAnnotationTable
%-------------------------------------------------------------------------------
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 3
    LDthreshold = 0.5;
end
%-------------------------------------------------------------------------------

numGenesChar = length(genesChar);
%-------------------------------------------------------------------------------
% Initialize:
numGWASMapped = zeros(numGenesChar,1);
numLDSNPs = zeros(numGenesChar,1);

for i = 1:numGenesChar
    gene_i = genesChar{i};

    %---------------------------------------------------------------------------
    % --numGWASMapped: the number of GWAS SNPs mapped directly to the gene
    %---------------------------------------------------------------------------
    isMappedGene = cellfun(@(x)f_isin(gene_i,x),SNPAnnotationTable.mappedGenes);
    theMappedSNPs = unique(SNPAnnotationTable.SNP(isMappedGene));
    numGWASMapped(i) = length(theMappedSNPs);

    %---------------------------------------------------------------------------
    % --numLD: the number of GWAS SNPs LD to a gene
    %---------------------------------------------------------------------------
    % (i) find SNPs for this gene:
    SNP_list = SQL_SNPsForGene(gene_i);
    numSNPs = length(SNP_list);

    % (ii) Get a list of all unique SNPs that are LD to these SNPs:
    LD_SNPs = cell(numSNPs,1);
    for j = 1:numSNPs
        mySNP = SNP_list{j};
        LD_SNPs{j} = SQL_SNP_LD_SNP(mySNP,LDthreshold);
    end
    all_LD_SNPs = unique(vertcat(LD_SNPs{:}));

    % (iii) Count how many match the list of disease-related SNPs:
    numLDSNPs = sum(ismember(all_LD_SNPs,SNPAnnotationTable.SNP));

    % If I want to count the number of SNPs that are LD to a given gene, I need to:
    % (i) Check every SNP in that gene (?? Using SNP_identifier.csv??)
    % (ii) Get a list of all (unique) SNPs are LD to these SNPs (Using LD_r2.csv)
    % (iii) Check all LD SNPs for direct disease annotations (in SNP_identifier.csv).
    % isGeneAndLD = strcmp(SNPAnnotationTable.mappedGene,gene_i) & SNPAnnotationTable.isLD;
    % theLDSNPs = unique(SNPAnnotationTable.SNP_id(isGeneAndLD));
    % NB: No way to double check that all LD annotations are from unique SNPs...
    % (is an assumption that Janette has formatted the file in this way)

    % Instead, we count all unique SNPs that are in a gene and LD to some GWAS hit:
    % if isempty(theLDSNPs)
    %     numLDSNPs(i) = 0;
    % else
    %     % Exclude SNPs that are mapped directly
    %     inBoth = intersect(theMappedSNPs,theLDSNPs);
    %     numLDSNPs(i) = length(theLDSNPs) - length(inBoth);
    % end
end

    function out = f_isin(A,B)
        if iscell(B)
            out = ismember(A,B);
        else
            out = strcmp(A,B);
        end
    end
end
