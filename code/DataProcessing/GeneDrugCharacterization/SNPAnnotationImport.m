function SNPAnnotationTable = SNPAnnotationImport(whatDisease,LDthreshold)
% Import data on SNP-gene relationships, and disease, GWAS/LD annotations
%-------------------------------------------------------------------------------

if nargin < 1
    whatDisease = 'all'; % don't filter
end
if nargin < 2
    LDthreshold = 0.5;
end

%-------------------------------------------------------------------------------
% Load processed SNP annotations generated using SNPAnnotationProcessing
fileName = sprintf('SNPAnnotationTable_%u.mat',LDthreshold*100);
load(fileName,'SNPAnnotationTable');

%-------------------------------------------------------------------------------
% Filter by disease:
%-------------------------------------------------------------------------------
sizePreFilter = height(SNPAnnotationTable);
switch whatDisease
case 'SZP'
    SNPAnnotationTable = SNPAnnotationTable(SNPAnnotationTable.isSZP,:);
case 'ADHD'
    SNPAnnotationTable = SNPAnnotationTable(SNPAnnotationTable.isADHD,:);
case 'ASD'
    SNPAnnotationTable = SNPAnnotationTable(SNPAnnotationTable.isASD,:);
case 'BIP'
    SNPAnnotationTable = SNPAnnotationTable(SNPAnnotationTable.isBIP,:);
case 'MDD'
    SNPAnnotationTable = SNPAnnotationTable(SNPAnnotationTable.isMDD,:);
case 'diabetes'
    SNPAnnotationTable = SNPAnnotationTable(SNPAnnotationTable.isDiabetes,:);
case 'all'
    % ---Keep all---
otherwise
    error('Unknown disease: ''%s''',whatDisease);
end
numAnnotations = height(SNPAnnotationTable);
if strcmp(whatDisease,'all')
    fprintf(1,'Keeping annotations from all diseases!\n');
else
    fprintf(1,'Filtering on %s reduced from %u -> %u annotations\n',...
                            whatDisease,sizePreFilter,numAnnotations);
end



% %-------------------------------------------------------------------------------
% % Match SNPs to mapped and LD-mapped genes using precomputed matching:
% %-------------------------------------------------------------------------------
% fileName = sprintf('SNP_gene_map_%u.mat',LDthreshold*100);
% load(fileName,'SNPgeneTable');
% % Match:
% [~,ia,ib] = intersect(SNPAnnotationTable.SNP_id,SNPgeneTable.uniqueSNPs);
% SNPAnnotationTable = SNPAnnotationTable(ia,:);
% SNPAnnotationTable.mappedGenes = SNPgeneTable.mappedGenes(ib);
% SNPAnnotationTable.LDgenes = SNPgeneTable.LDgenes(ib);
%
% %-------------------------------------------------------------------------------
% % Give outputs to screen:
% %-------------------------------------------------------------------------------
% fprintf(1,'%u annotations for %u GWAS mapped and %u LD\n',numAnnotations,...
%                     sum(SNPAnnotationTable.isGWAS),sum(SNPAnnotationTable.isLD));
%
% hasGeneName = cellfun(@(x)~isempty(x),SNPAnnotationTable.mappedGene);
% fprintf(1,'%u/%u annotations are mapped to gene names\n',sum(hasGeneName),numAnnotations);
% fprintf(1,'%u/%u GWAS-mapped annotations have gene names\n',...
%                 sum(SNPAnnotationTable.isGWAS & hasGeneName),sum(SNPAnnotationTable.isGWAS));
% fprintf(1,'%u/%u LD annotations have gene names\n',...
%             sum(SNPAnnotationTable.isLD & hasGeneName),sum(SNPAnnotationTable.isLD));
%
% % All SNPs with a disease annotation (GWAS or LD):
% allDiseaseSNPs = unique(SNPAnnotationTable.SNP_id);
%
% %-------------------------------------------------------------------------------
% % Generate a SNP -> Gene map
% %-------------------------------------------------------------------------------
% SNPGeneMap = SNPAnnotationTable(cellfun(@(x)~isempty(x),SNPAnnotationTable.mappedGene),1:2);
% % Make unique:
% [~,ix] = unique(SNPGeneMap.SNP_id);
% SNPGeneMap = SNPGeneMap(ix,:);
% fprintf(1,'%u/%u SNPs have a matching gene name (%u genes)\n',height(SNPGeneMap),...
%         length(unique(SNPAnnotationTable.SNP_id)),length(unique(SNPGeneMap.mappedGene)));

end
