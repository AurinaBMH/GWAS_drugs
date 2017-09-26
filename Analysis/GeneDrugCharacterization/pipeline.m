% Pipeline for producing table characterizing individual genes

%-------------------------------------------------------------------------------
% Get genes
%-------------------------------------------------------------------------------
% List of all genes with drug targets in Drugbank are in **gene_ATC_matrix.csv**

fid = fopen('8_1_drug_gene.csv','r');
C = textscan(fid,'%s%s','Delimiter',',');
fclose(fid)
geneName = C{1}; drugName = C{2};
geneDrug = table(geneName,drugName);
allUniqueGenes = unique(geneDrug.geneName);
numUniqueGenes = length(allUniqueGenes);
fprintf(1,'%u unique genes with drug targets (in Drugbank)\n',numUniqueGenes);
numUniqueDrugs = length(unique(geneDrug.drugName));
fprintf(1,'%u unique drugs with gene targets (in Drugbank)\n',numUniqueDrugs);



%-------------------------------------------------------------------------------
