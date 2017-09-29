function [eQTLTable,isEGene,isSNPGene] = eQTLImport()
% Import data on eQTL-SNP relationships
%-------------------------------------------------------------------------------

fid = fopen('4_GWAS_eQTL_combn.csv','r');
C = textscan(fid,'%s%s%f%s%f','Delimiter',',','HeaderLines',1);
fclose(fid);
eQTL = C{1};
SNP = C{2};
eGene_qVal = C{3};
eGene_qVal(eGene_qVal==-999) = NaN;
tissueSample = C{4};
SNPgene_pVal = C{5};
SNPgene_pVal(SNPgene_pVal==-999) = NaN;
eQTLTable = table(eQTL,SNP,eGene_qVal,tissueSample,SNPgene_pVal);
isMissing = strcmp(eQTL,'0');
fprintf(1,'%u eQTLs assigned ''0'' -- no hgnc gene symbol -- removed\n',sum(isMissing));
eQTLTable = eQTLTable(~isMissing,:);
fprintf(1,'%u eQTL annotations loaded\n',height(eQTLTable));
isEGene = ~isnan(eQTLTable.eGene_qVal);
fprintf(1,'%u eGene annotations\n',sum(isEGene));
isSNPGene = ~isnan(eQTLTable.SNPgene_pVal);
fprintf(1,'%u SNPgene annotations\n',sum(isSNPGene));

end
