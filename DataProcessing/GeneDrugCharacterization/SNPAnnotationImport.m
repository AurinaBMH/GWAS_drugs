function [SNPAnnotationTable,SNPGeneMap,allDiseaseSNPs] = SNPAnnotationImport()
% Import data on SNP-gene relationships, and disease, GWAS/LD annotations
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
mappedGene = C{7};
mappedGene(strcmp(mappedGene,'0')) = {''}; % remove '0' -> empty
isGWAS = logical(C{8});
isLD = logical(C{9});
SNPAnnotationTable = table(SNP_id,mappedGene,isGWAS,isLD,isSZP,isADHD,isASD,isBIP,isMDD);
%-------------------------------------------------------------------------------

numAnnotations = height(SNPAnnotationTable);
fprintf(1,'%u annotations for %u GWAS mapped and %u LD\n',numAnnotations,sum(isGWAS),sum(isLD));

hasGeneName = cellfun(@(x)~isempty(x),mappedGene);
fprintf(1,'%u/%u annotations are mapped to gene names\n',sum(hasGeneName),numAnnotations);
fprintf(1,'%u/%u GWAS-mapped annotations have gene names\n',sum(isGWAS & hasGeneName),sum(isGWAS));
fprintf(1,'%u/%u LD annotations have gene names\n',sum(isLD & hasGeneName),sum(isLD));

% Cell of all SNPs with a disease annotation (GWAS or LD):
allDiseaseSNPs = unique(SNPAnnotationTable.SNP_id);

%-------------------------------------------------------------------------------
% Generate a SNP->Gene map
%-------------------------------------------------------------------------------
SNPGeneMap = SNPAnnotationTable(cellfun(@(x)~isempty(x),SNPAnnotationTable.mappedGene),1:2);
% Make unique:
[~,ix] = unique(SNPGeneMap.SNP_id);
SNPGeneMap = SNPGeneMap(ix,:);
fprintf(1,'%u/%u SNPs have a matching gene name (%u genes)\n',height(SNPGeneMap),...
        length(unique(SNPAnnotationTable.SNP_id)),length(unique(SNPGeneMap.mappedGene)));

end
