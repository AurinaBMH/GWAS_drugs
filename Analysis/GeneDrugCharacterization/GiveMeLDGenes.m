function theLDGenes = GiveMeLDGenes(myGene,SNPGeneMap,LDRelateTable,allDiseaseSNPs)
%-------------------------------------------------------------------------------
% 1. Find SNPs for a given gene (from SNPAnnotationTable)
% 2. Find SNPs LD to those SNPs
% 3. Match LD SNPs to genes -> output
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% --1-- Get SNPs for the gene:
theSNPs = SNPGeneMap.SNP_id(strcmp(SNPGeneMap.mappedGene,myGene));
if isempty(theSNPs)
    theLDGenes = {};
    return
end

%-------------------------------------------------------------------------------
% --2-- Get SNPs LD to these SNPs (that have a relevant disease annotation):
numSNPs = length(theSNPs);
theLDSNPs = cell(numSNPs,1);
for i = 1:numSNPs
    theLDSNPs{i} = GiveMeLDSNPs(theSNPs{i},LDRelateTable,allDiseaseSNPs);
end
theLDSNPs = vertcat(theLDSNPs{:});
if isempty(theLDSNPs)
    theLDGenes = {};
    return
end

%-------------------------------------------------------------------------------
% --3-- Map SNPs to genes
theLDSNPs = unique(theLDSNPs);
theLDGenes = cellfun(@(x)SNPGeneMap.mappedGene(strcmp(SNPGeneMap.SNP_id,x)),theLDSNPs,'UniformOutput',false);
theLDGenes = unique(vertcat(theLDGenes{:}));

if ismember(myGene,theLDGenes)
    theLDGenes = theLDGenes(~strcmp(theLDGenes,myGene));
    fprintf(1,'Removed target gene, %s, from theLDGene list\n',myGene);
end

end
