function theLDGenes = GiveMeLDGenes(myGene,threshold); %SNPGeneMap,LDRelateTable,allDiseaseSNPs)
%-------------------------------------------------------------------------------
% 1. Find SNPs for a given gene (SQL_SNPsForGene)
% 2. Find SNPs LD to those SNPs (SQL_SNP_LD_SNP)
% 3. Match LD SNPs to genes -> output (SQL_genesForSNPs)
%-------------------------------------------------------------------------------

if nargin < 2
    threshold = 0.5;
end

theLDGenes = {};

% --1-- Get SNPs for the gene:
gene_SNPs = SQL_SNPsForGene(geneName);
if isempty(gene_SNPs), return; end

%-------------------------------------------------------------------------------
% --2-- Get all SNPs LD to these SNPs:
numSNPs = length(gene_SNPs);
theLDSNPs = cell(numSNPs,1);
for i = 1:numSNPs
    theLDSNPs{i} = SQL_SNP_LD_SNP(gene_SNPs{i},threshold);
end
theLDSNPs = unique([theLDSNPs{:}]);
if isempty(theLDSNPs), return; end

%-------------------------------------------------------------------------------
% --3-- Map SNPs to genes
theLDGenes = SQL_genesForSNPs(theLDSNPs);

% Remove the target gene from the list...?
if ismember(myGene,theLDGenes)
    theLDGenes = theLDGenes(~strcmp(theLDGenes,myGene));
    fprintf(1,'Removed target gene, %s, from theLDGene list\n',myGene);
end

% sqlMethod = true;
% if sqlMethod
%
% else
%     %---------------------------------------------------------------------------
%     % Have to use restricted data subset:
%     %---------------------------------------------------------------------------
%     % --1-- Get SNPs for the gene:
%     theSNPs = SNPGeneMap.SNP_id(strcmp(SNPGeneMap.mappedGene,myGene));
%     if isempty(theSNPs)
%         theLDGenes = {};
%         return
%     end
%
%     %-------------------------------------------------------------------------------
%     % --2-- Get SNPs LD to these SNPs (that have a relevant disease annotation):
%     numSNPs = length(theSNPs);
%     theLDSNPs = cell(numSNPs,1);
%     for i = 1:numSNPs
%         theLDSNPs{i} = GiveMeLDSNPs(theSNPs{i},LDRelateTable,allDiseaseSNPs);
%     end
%     theLDSNPs = vertcat(theLDSNPs{:});
%     if isempty(theLDSNPs)
%         theLDGenes = {};
%         return
%     end
%
%     %-------------------------------------------------------------------------------
%     % --3-- Map SNPs to genes
%     theLDSNPs = unique(theLDSNPs);
%     theLDGenes = cellfun(@(x)SNPGeneMap.mappedGene(strcmp(SNPGeneMap.SNP_id,x)),theLDSNPs,'UniformOutput',false);
%     theLDGenes = unique(vertcat(theLDGenes{:}));
%
%     if ismember(myGene,theLDGenes)
%         theLDGenes = theLDGenes(~strcmp(theLDGenes,myGene));
%         fprintf(1,'Removed target gene, %s, from theLDGene list\n',myGene);
%     end
% end

end
