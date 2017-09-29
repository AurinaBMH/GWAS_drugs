function theLDSNPs = GiveMeLDSNPs(mySNP,LDRelateTable,allDiseaseSNPs)
%-------------------------------------------------------------------------------
% Counts how many unique SNPs are LD to a given SNP...
%-------------------------------------------------------------------------------

% Get all SNPs LD to the input SNP:
theLDSNP1 = LDRelateTable.SNP_id_2(strcmp(LDRelateTable.SNP_id_1,mySNP));
theLDSNP2 = LDRelateTable.SNP_id_1(strcmp(LDRelateTable.SNP_id_2,mySNP));
theLDSNPs = unique(vertcat(theLDSNP1,theLDSNP2));

if isempty(theLDSNPs)
    return
end

% Check which are in the SNPAnnotationTable:
isInAnnotation = ismember(theLDSNPs,allDiseaseSNPs);
if ~all(isInAnnotation)
    warning('Not all LD SNPs have annotations?')
end

% Count how many:
theLDSNPs = theLDSNPs(isInAnnotation);

end
