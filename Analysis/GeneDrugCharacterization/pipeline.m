% Pipeline for producing table characterizing individual genes

%===============================================================================
% LOAD AND PROCESS DATA
%===============================================================================

%-------------------------------------------------------------------------------
% Get genes
%-------------------------------------------------------------------------------
% List of all genes with drug targets in Drugbank are in **gene_ATC_matrix.csv**

fid = fopen('8_1_drug_gene.csv','r');
C = textscan(fid,'%s%s','Delimiter',',');
fclose(fid);
geneName = C{1}; drugName = C{2};
geneDrug = table(geneName,drugName);
allUniqueGenes = unique(geneDrug.geneName);
numUniqueGenes = length(allUniqueGenes);
fprintf(1,'%u unique genes with drug targets (in Drugbank)\n',numUniqueGenes);
numUniqueDrugs = length(unique(geneDrug.drugName));
fprintf(1,'%u unique drugs with gene targets (in Drugbank)\n',numUniqueDrugs);

%-------------------------------------------------------------------------------
% Read in SNP annotations
%-------------------------------------------------------------------------------
fid = fopen('2_1_SNP_identifier.csv','r');
C = textscan(fid,'%s%u%u%u%u%u%s%u%u','Delimiter',',','HeaderLines',1);
fclose(fid);
SNP_id = C{1};
isSZP = logical(C{2});
isADHD = logical(C{3});
isASD = logical(C{4});
isBIP = logical(C{5});
isMDD = logical(C{6});
mappedGeneSNP = C{7};
mappedGeneSNP(strcmp(mappedGeneSNP,'0')) = {''}; % remove '0' -> empty
isGWAS = logical(C{8});
isLD = logical(C{9});
SNPAnnotationTable = table(SNP_id,mappedGeneSNP,isGWAS,isLD,isSZP,isADHD,isASD,isBIP,isMDD);

%-------------------------------------------------------------------------------
% Now loop through genes to characterize...
numGWASMapped = zeros(numUniqueGenes,1);

for i = 1:numUniqueGenes
    fprintf(1,'[%u/%u]: %s\n',i,numUniqueGenes,allUniqueGenes{i});
    % --numGWASMapped: the number of GWAS hits mapped directly to a gene

end
