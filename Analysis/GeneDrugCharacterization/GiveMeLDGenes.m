function theLDgenes = GiveMeLDGenes(myGene,SNPGeneMap,LDRelateTable,allDiseaseSNPs)
%-------------------------------------------------------------------------------
% 1. Find SNPs for a given gene (from SNPAnnotationTable)
% 2. Find SNPs LD to those SNPs
% 3. Match LD SNPs to genes -> output
%-------------------------------------------------------------------------------

% --1-- Get SNPs for the gene:
theSNPs = SNPGeneMap.SNP_id(strcmp(SNPGeneMap.mappedGene,myGene));

if isempty(theSNPs)
    theLDGenes = {};
    return
end

% --2-- Get SNPs LD to these SNPs:
numSNPs = length(theSNPs);
theLDSNPs = cell(numSNPs,1);
for i = 1:numSNPs
    theLDSNPs{i} = GiveMeLDSNPs(theSNPs{i},LDRelateTable,allDiseaseSNPs)
end


theLDGenes = cellfun(@(x)SNPGeneMap.mappedGene(strcmp(SNPGeneMap.SNP_id,x)),theLDSNPs,'UniformOutput',false);
theSNPs =

% Get all SNPs LD to the input SNP:
theLDSNP1 = LDRelateTable.SNP_id_2(strcmp(LDRelateTable.SNP_id_1,mySNP));
theLDSNP2 = LDRelateTable.SNP_id_1(strcmp(LDRelateTable.SNP_id_2,mySNP));
theLDSNPs = unique(vertcat(theLDSNP1,theLDSNP2));

% if isempty(theLDSNPs)
%     numLD = 0;
%     return
% end

% Check which are in the SNPAnnotationTable:
isInAnnotation = ismember(theLDSNPs,SNPAnnotationTable.SNP_id);

if ~all(isInAnnotation)
    warning('Not all LD SNPs have annotations?')
end

% Count how many:
theLDSNPs = theLDSNPs(isInAnnotation);

end
