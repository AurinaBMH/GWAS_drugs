function DNAStats = TellMeDNADistance(genesChar,SNPAnnotationTable,LDthreshold)
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
numSNPs = zeros(numGenesChar,1);
numGWAS = zeros(numGenesChar,1);
numLD = zeros(numGenesChar,1);
percGWAS = zeros(numGenesChar,1);
percLD = zeros(numGenesChar,1);

fprintf(1,'Computing number of SNPs (mapped and LD) to SNPs in %u genes\n',numGenesChar);
for i = 1:numGenesChar
    gene_i = genesChar{i};
    fprintf(1,'%u/%u: %s',i,numGenesChar,gene_i);

    %---------------------------------------------------------------------------
    % --numSNPs: the number of SNPs annotated to the gene (how big it is?)
    %---------------------------------------------------------------------------
    SNP_list = SQL_SNPsForGene(gene_i);
    numSNPs(i) = length(SNP_list);

    %---------------------------------------------------------------------------
    % --numGWASMapped: the number of GWAS SNPs mapped directly to the gene
    %---------------------------------------------------------------------------
    isMappedGene = cellfun(@(x)ismember(gene_i,x),SNPAnnotationTable.mappedGenes);
    theMappedSNPs = unique(SNPAnnotationTable.SNP(isMappedGene));
    numGWAS(i) = length(theMappedSNPs);

    %---------------------------------------------------------------------------
    % --numLD: the number of GWAS SNPs LD to this gene
    %---------------------------------------------------------------------------
    isLDGene = cellfun(@(x)ismember(gene_i,x),SNPAnnotationTable.LDgenes);
    theLDSNPs = unique(SNPAnnotationTable.SNP(isLDGene));
    numLD(i) = length(theLDSNPs);

    %---------------------------------------------------------------------------
    % Normalize:
    %---------------------------------------------------------------------------
    percGWAS(i) = numGWAS(i)/numSNPs(i);
    percLD(i) = numLD(i)/numSNPs(i);

    %---------------------------------------------------------------------------
    % --numLD_2: the number of GWAS SNPs LD to a gene
    %---------------------------------------------------------------------------
    % % (i) find SNPs for this gene:
    % SNP_list = SQL_SNPsForGene(gene_i);
    % numSNPs = length(SNP_list);
    %
    % % (ii) Get a list of all unique SNPs that are LD to these SNPs:
    % LD_SNPs = cell(numSNPs,1);
    % for j = 1:numSNPs
    %     mySNP = SNP_list{j};
    %     LD_SNPs{j} = SQL_SNP_LD_SNP(mySNP,LDthreshold);
    % end
    % all_LD_SNPs = unique(vertcat(LD_SNPs{:}));
    %
    % isMappedSNP = ismember(all_LD_SNPs,theMappedSNPs);
    % all_LD_SNPs = all_LD_SNPs(~isMappedSNP);
    %
    % % (iii) Count how many match the list of disease-related SNPs:
    % % Probably want to exclude SNPs used above
    % numLDSNPs(i) = sum(ismember(all_LD_SNPs,SNPAnnotationTable.SNP));
    % %
    fprintf(1,' (%u,%u / %u)\n',numGWAS(i),numLD(i),numSNPs(i));

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

%-------------------------------------------------------------------------------
% Output structure
DNAStats = struct();
DNAStats.numSNPs = numSNPs;
DNAStats.numGWAS = numGWAS;
DNAStats.percGWAS = percGWAS;
DNAStats.numLD = numLD;
DNAStats.percLD = percLD;

    % function out = f_isin(A,B)
    %     if iscell(B)
    %         out = ismember(A,B);
    %     else
    %         out = strcmp(A,B);
    %     end
    % end
end
