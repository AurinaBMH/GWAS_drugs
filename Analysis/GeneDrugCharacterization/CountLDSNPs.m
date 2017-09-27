function numLDSNPs = CountLDSNPs(theGene,SNPAnnotationTable,LDRelateTable,SNPGeneMap,beVocal)
%-------------------------------------------------------------------------------
% Counts the number of SNPs LD to a given target SNP with disease annotations
%-------------------------------------------------------------------------------
if nargin < 5
    beVocal = false;
end

% (i) Check every SNP in the gene (Maybe all are in SNPAnnotationTable???)
allSNPs = SNPGeneMap.SNP_id(strcmp(SNPGeneMap.mappedGeneSNP,theGene));
numSNPs = length(allSNPs);

if beVocal
    fprintf(1,'%u SNPs for %s\n',numSNPs,theGene);
end

if numSNPs==0
    % No SNPs in this gene?? So cannot be any LD matches:
    numLDSNPs = 0;
    return
end

% (ii) Get a list of all (unique) SNPs are LD to these SNPs (Using LDRelateTable)
theLDSNPs = cell(numSNPs,1);
for j = 1:numSNPs
    theLDSNP1 = LDRelateTable.SNP_id_2(strcmp(LDRelateTable.SNP_id_1,SNPAnnotationTable.SNP_id{j}));
    theLDSNP2 = LDRelateTable.SNP_id_1(strcmp(LDRelateTable.SNP_id_2,SNPAnnotationTable.SNP_id{j}));
    theLDSNPs{j} = unique(vertcat(theLDSNP1,theLDSNP2));
end
% Now agglomerate:
allLDSNPs = unique(vertcat(theLDSNPs{:}));
numLDSNPs = length(allLDSNPs);

% (iii) Check all LD SNPs for direct disease annotations (in SNP_identifier.csv).


end
