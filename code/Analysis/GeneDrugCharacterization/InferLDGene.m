%-------------------------------------------------------------------------------
% Actually this was a bit of a failure... :-/
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
% Infer the LD gene (i.e., the gene causing the annotation) for LD annotations
% in SNPAnnotationTable
%-------------------------------------------------------------------------------
% For LD SNPs, we can annotate the SNP that caused the annotation
% LDGene = cell(numAnnotations,1);
% for i = 1:numAnnotations
%     if SNPAnnotationTable.isLD(i)
%         % We want to determine which gene caused the annotation:
%         theLDSNP1 = LDRelateTable.SNP_id_2(strcmp(LDRelateTable.SNP_id_1,SNPAnnotationTable.SNP_id{i}));
%         theLDSNP2 = LDRelateTable.SNP_id_1(strcmp(LDRelateTable.SNP_id_2,SNPAnnotationTable.SNP_id{i}));
%         theLDSNPs = unique(vertcat(theLDSNP1,theLDSNP2));
%         if ischar(theLDSNPs), theLDSNPs = {theLDSNPs}; end
%
%         % Match to genes:
%         theLDGenes = cellfun(@(x)SNPGeneMap.mappedGene(strcmp(SNPGeneMap.SNP_id,x)),theLDSNPs,'UniformOutput',false);
%
%         if isempty(theLDGenes)
%             % No matches? Must be an LD relationship to a SNP without a gene
%             warning('No LD gene for %s',SNPAnnotationTable.SNP_id{i});
%             LDGene{i} = '';
%         elseif length(theLDGenes)==1
%             % One match -- easy:
%             LDGene{i} = theLDGenes{1};
%         else
%             fprintf(1,'Multiple matches :-O\n');
%         end
%
%     else
%         LDGene{i} = '';
%     end
% end
