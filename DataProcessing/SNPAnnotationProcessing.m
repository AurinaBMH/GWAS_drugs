% function [SNPAnnotationTable,SNPgeneTable] = SNPAnnotationProcessing(LDthreshold)
% Idea is to make a nicer version of the SNPAnnotationTable

% if nargin < 1
    LDthreshold = 0.5;
% end

%-------------------------------------------------------------------------------
% IMPORT GWAS-BASED INFORMATION:
%-------------------------------------------------------------------------------
fid = fopen('2_1_SNP_identifier_v2.csv','r');
fprintf(1,'Using new file generated Jul-2018\n')
C = textscan(fid,'%s%u%u%u%u%u%s%u%u%u','Delimiter',',','HeaderLines',1);
fclose(fid);
SNP_id = C{1};
isSZP = logical(C{2});
isADHD = logical(C{3});
isASD = logical(C{4});
isBIP = logical(C{5});
isMDD = logical(C{6});
isDiabetes = logical(C{10});
mappedGene = C{7};
mappedGene(strcmp(mappedGene,'0')) = {''}; % remove '0' -> empty
isGWAS = logical(C{8});
isLD = logical(C{9});
% SNPAnnotationTable = table(SNP_id,mappedGene,isGWAS,isLD,isSZP,isADHD,isASD,isBIP,isMDD,isDiabetes);
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% We first want a table relating each SNP to a gene
% ONLY INCLUDING GWAS:
%-------------------------------------------------------------------------------
mappedGeneJanette = mappedGene;
SNPAnnotationTableAll = table(SNP_id,isSZP,isADHD,isASD,isBIP,...
                                    isMDD,isDiabetes,mappedGeneJanette);

%-------------------------------------------------------------------------------
% Filter to list only GWAS-annotated SNPs:
SNPAnnotationTableAll = SNPAnnotationTableAll(isGWAS,:);

%-------------------------------------------------------------------------------
% Filter to unique SNPs--—SNPAnnotationTable
[SNP,ia] = unique(SNPAnnotationTableAll.SNP_id);
numUniqueSNPs = length(SNP);
fprintf(1,'%u unique GWAS SNPs have been implicated in disease\n',numUniqueSNPs);
isSZP = false(numUniqueSNPs,1);
isADHD = false(numUniqueSNPs,1);
isASD = false(numUniqueSNPs,1);
isBIP = false(numUniqueSNPs,1);
isMDD = false(numUniqueSNPs,1);
isDiabetes = false(numUniqueSNPs,1);
for i = 1:numUniqueSNPs
    theSNP = SNP{i};
    tableSNP = SNPAnnotationTableAll(strcmp(SNPAnnotationTableAll.SNP_id,theSNP),:);
    if any(tableSNP.isSZP)
        isSZP(i) = true;
    end
    if any(tableSNP.isADHD)
        isADHD(i) = true;
    end
    if any(tableSNP.isASD)
        isASD(i) = true;
    end
    if any(tableSNP.isBIP)
        isBIP(i) = true;
    end
    if any(tableSNP.isMDD)
        isMDD(i) = true;
    end
    if any(tableSNP.isDiabetes)
        isDiabetes(i) = true;
    end
end
SNPAnnotationTable = table(SNP,isSZP,isADHD,isASD,isBIP,isMDD,isDiabetes);

%-------------------------------------------------------------------------------
% For all unique SNPs, we want a column: LDgenes
%-------------------------------------------------------------------------------
% Initialize:
mappedGenes = cell(numUniqueSNPs,1);
LD_SNPs = cell(numUniqueSNPs,1);
LDgenes = cell(numUniqueSNPs,1);
for i = 1:numUniqueSNPs
    theSNP = SNP{i};
    fprintf(1,'%u/%u: %s\n',i,numUniqueSNPs,theSNP);

    % Get the mapped gene for this SNP:
    mappedGenes{i} = SQL_genesForSNPs(theSNP);

    % Get LD SNPs:
    LD_SNPs{i} = SQL_SNP_LD_SNP(theSNP,LDthreshold);

    % Map each SNP to its set of LD genes (if there are any to match):
    if ~isempty(LD_SNPs{i})
        LDgenes{i} = SQL_genesForSNPs(LD_SNPs{i});
    end
end
SNPgeneTable = table(SNP,mappedGenes,LD_SNPs,LDgenes);

%-------------------------------------------------------------------------------
% Save to a .mat file
fileName = fullfile('DataOutput',sprintf('SNP_gene_map_%u.mat',LDthreshold*100));
save(fileName,'SNPgeneTable','LDthreshold');
fprintf(1,'Saved SNP/LD/Gene mapping to %s\n',fileName);

%-------------------------------------------------------------------------------
% Match and integrate:
SNPAnnotationTable.mappedGenes = mappedGenes;
SNPAnnotationTable.LDgenes = LDgenes;

% noMatches = cellfun(@isempty,SNPAnnotationTable.mappedGenes) & cellfun(@isempty,SNPAnnotationTable.LDgenes);
% SNPAnnotationTable = SNPAnnotationTable(~noMatches,:);

% numEntries = height(SNPAnnotationTable);
% mappedGene = cell(numEntries,1);
% LDgenes = cell(numEntries,1);
% for i = 1:numEntries
%     theSNP = SNPAnnotationTable.SNP{i};
%     ind = find(strcmp(SNPgeneTable.SNP,theSNP));
%     if ~isempty(ind)
%         mappedGene{i} = SNPgeneTable.mappedGenes{ind};
%         LDgenes{i} = SNPgeneTable.LDgenes{ind};
%     end
% end
% SNPAnnotationTable.mappedGene = mappedGene;
% SNPAnnotationTable.LDgenes = LDgenes;

%-------------------------------------------------------------------------------
% Save to a .mat file
fileName = fullfile('DataOutput',sprintf('SNPAnnotationTable_%u.mat',LDthreshold*100));
save(fileName,'SNPAnnotationTable');
fprintf(1,'Saved SNP annotation table to %s\n',fileName);

% end
