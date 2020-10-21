% check AHBA measure
clear all; 
load('GWAS_disordersMAGMA.mat')
whatDisease = 'SCZ'; 

params = SetDefaultParams();
whatDrugTargets = params.whatDrugTargets; 

indicatorTable = ImportTreatmentLists(false, whatDrugTargets);
allTargetGenes = indicatorTable.Properties.RowNames;
numTargetGenes = length(allTargetGenes);
fprintf(1,'Analyzing %u genes that have drug targets in our list\n',numTargetGenes);

%-------------------------------------------------------------------------------
% Get genes for a given GWAS study for different mapping methods:
%-------------------------------------------------------------------------------

listGENESmapped = DISORDERlist.MAGMAdefault.(whatDisease);

pThr_m = 0.05/size(listGENESmapped,1); % Bonf correction for the number of genes in the list     
allGWASgenes = listGENESmapped.GENENAME(listGENESmapped.P<pThr_m);

% load gene expression data
[geneCoexp,AllenGeneInfo] = LoadCoexpression();
isContext = ismember(AllenGeneInfo.GeneSymbol,allGWASgenes);
AllenMeanCoexp = nan(numTargetGenes,1);

for i = 1:numTargetGenes
    
    gene_i = allTargetGenes{i};

    allenIndex = strcmp(AllenGeneInfo.GeneSymbol,gene_i);
    if any(allenIndex)
        % Compute the coexpression values of this gene to the set of context genes:
        % take the absolute values of coexpression
        V = geneCoexp(allenIndex,isContext); 
%         figure; 
%         subplot(1,2,1); histogram(V, 10); 
%         title(sprintf('%s', gene_i))
%         xlabel(sprintf('%s and GWAS genes coexpression', gene_i))
%         
%         subplot(1,2,2); histogram(abs(V), 10); 
%         title(sprintf('%s', gene_i))
%         xlabel(sprintf('%s and GWAS genes abs coexpression', gene_i))
        
        coExpContext = abs(V);
        AllenMeanCoexp(i) = nanmean(coExpContext);
    else
        % This gene could not be matched to AHBA data
        % warning('%s could not be matched to the Allen expression data',gene_i)
        % AllenMeanCoexp(i) = NaN;
    end
end

% if there are no real differences, then by selecting GWAS genes for each
% disorder, we select random sets of genes - it's expected that overall,
% random sets of genes will give similar coexpression values. 
figure; 
subplot(1,2,1); histogram(V_diabetes); 
subplot(1,2,2); histogram(V_scz); 



geneScores.AllenMeanCoexpMapped = TellMeAllenCoexp(allTargetGenes,allGWASgenes);