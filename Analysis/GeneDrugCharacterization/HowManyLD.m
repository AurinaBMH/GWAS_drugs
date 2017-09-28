function numLD = HowManyLD(mySNP,SNPAnnotationTable,LDRelateTable)
%-------------------------------------------------------------------------------
% Counts how many unique SNPs are LD to a given SNP...
%-------------------------------------------------------------------------------

% Get all LD SNPs:
theLDSNP1 = LDRelateTable.SNP_id_2(strcmp(LDRelateTable.SNP_id_1,mySNP));
theLDSNP2 = LDRelateTable.SNP_id_1(strcmp(LDRelateTable.SNP_id_2,mySNP));
theLDSNPs = unique(vertcat(theLDSNP1,theLDSNP2));

if isempty(theLDSNPs)
    numLD = 0;
    return
end

% Check which are in the SNPAnnotationTable:
isInAnnotation = ismember(theLDSNPs,SNPAnnotationTable.SNP_id);

% Count how many:
numLD = sum(isInAnnotation);

end
