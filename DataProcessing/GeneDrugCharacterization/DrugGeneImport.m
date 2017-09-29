function [geneDrugTable,allUniqueGenes] = DrugGeneImport()
% Import drug-gene relationships
%-------------------------------------------------------------------------------

fid = fopen('8_1_drug_gene.csv','r');
C = textscan(fid,'%s%s','Delimiter',',');
fclose(fid);
geneName = C{1}; drugName = C{2};
geneDrugTable = table(geneName,drugName);
allUniqueGenes = unique(geneDrug.geneName);
numUniqueGenes = length(allUniqueGenes);
fprintf(1,'%u genes have drug targets in Drugbank\n',numUniqueGenes);
numUniqueDrugs = length(unique(geneDrug.drugName));
fprintf(1,'%u drugs have gene targets in Drugbank\n',numUniqueDrugs);

end
