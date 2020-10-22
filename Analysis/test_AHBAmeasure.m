% check AHBA measure
clear all;close all; 

whatDiseases = {'ADHD', 'BIP2', 'MDD2','SCZ', 'DIABETES'};
load('GWAS_disordersMAGMA.mat')
numGenes = 5; 

params = SetDefaultParams();

indicatorTable = ImportTreatmentLists(false, params.whatDrugTargets);
allTargetGenes = indicatorTable.Properties.RowNames;
numTargetGenes = length(allTargetGenes);
fprintf(1,'Analyzing %u genes that have drug targets in our list\n',numTargetGenes);

% load gene expression data
[geneCoexp,AllenGeneInfo] = LoadCoexpression();

% plot coexpression of several target genes to all other genes - are they
% different? 
figure('color','w', 'Position', [300, 300, 1500, 600]);
for g = 1:numGenes
    
   gene_g = allTargetGenes{g}; 
   allenIndex = strcmp(AllenGeneInfo.GeneSymbol,gene_g);
   
   if any(allenIndex)
       subplot(2,numGenes,g); 
       histogram(geneCoexp(allenIndex,:))
       title(sprintf('%s', allTargetGenes{g}))
       xlabel('coexpression')
       xlim([-1 1])
        
        
       subplot(2,numGenes,g+numGenes); 
       histogram(abs(geneCoexp(allenIndex,:)))
       title(sprintf('%s', allTargetGenes{g}))
       xlabel('abs(coexpression)')
       xlim([0 1])
   end
end


% now visualise coexpression to GWAS-selected genes for each disorder


for g = 1:numGenes
    % separate plot for each gene
    figure('color','w', 'Position', [300, 300, 1500, 600]);
    
    gene_g = allTargetGenes{g};
    allenIndex = strcmp(AllenGeneInfo.GeneSymbol,gene_g);
    
    for d=1:length(whatDiseases)
        listGENESmapped = DISORDERlist.MAGMAdefault.(whatDiseases{d});
        
        pThr_m = 0.05/size(listGENESmapped,1); % Bonf correction for the number of genes in the list
        allGWASgenes = listGENESmapped.GENENAME(listGENESmapped.P<pThr_m);
        isGWAS = ismember(AllenGeneInfo.GeneSymbol,allGWASgenes);
        
        if any(allenIndex)
            subplot(2,length(whatDiseases),d);
            histogram(geneCoexp(allenIndex,isGWAS), 10)
            title(sprintf('GWAS: %s, gene %s', whatDiseases{d}, allTargetGenes{g}))
            xlabel(sprintf('coexpression, mean = %.4f', nanmean(geneCoexp(allenIndex,isGWAS))))
            xlim([-1 1])
            
            subplot(2,length(whatDiseases),d+length(whatDiseases));
            histogram(abs(geneCoexp(allenIndex,isGWAS)), 10)
            title(sprintf('GWAS: %s, gene %s', whatDiseases{d}, allTargetGenes{g}))
            xlabel(sprintf('abs(coexpression), mean = %.4f', nanmean(abs(geneCoexp(allenIndex,isGWAS)))))
            xlim([0 1])
        end
        
    end
end


    
%-------------------------------------------------------------------------------
% Get genes for a given GWAS study for different mapping methods:
%-------------------------------------------------------------------------------
AllenMeanCoexp = nan(length(whatDiseases), numTargetGenes);
AllenMeanCoexp_abs = nan(length(whatDiseases), numTargetGenes);

for d=1:length(whatDiseases)
    
    listGENESmapped = DISORDERlist.MAGMAdefault.(whatDiseases{d});
    pThr_m = 0.05/size(listGENESmapped,1); % Bonf correction for the number of genes in the list
    allGWASgenes = listGENESmapped.GENENAME(listGENESmapped.P<pThr_m);
    isGWAS = ismember(AllenGeneInfo.GeneSymbol,allGWASgenes);
    
    
    for i = 1:numTargetGenes
        
        gene_i = allTargetGenes{i};
        allenIndex = strcmp(AllenGeneInfo.GeneSymbol,gene_i);
        
        if any(allenIndex)
            % Compute the coexpression values of this gene to the set of context genes:
            % take the absolute values of coexpression
            V = geneCoexp(allenIndex,isGWAS);

            AllenMeanCoexp_abs(d,i) = nanmean(abs(V));
            AllenMeanCoexp(d,i) = nanmean(V);
        end
    end
end

% how correlated absolute and regular scores are across disorders? 
for d1=1:length(whatDiseases)
    for d2=d1+1:length(whatDiseases)
        
        r_abs(d1, d2) = corr(AllenMeanCoexp_abs(d1,:)', AllenMeanCoexp_abs(d2,:)', 'rows', 'complete');
        r(d1, d2) = corr(AllenMeanCoexp(d1,:)', AllenMeanCoexp(d2,:)', 'rows', 'complete'); 
        
    end
end
[COL]=cbrewer('div', 'RdBu', 100);
colors = flipud(COL);

figure('color', 'w'); 
subplot(1,2,1); imagesc(r_abs); colormap(colors); caxis([-1 1]); colorbar; 
title('Correlation between abs(coexp)'); axis('square');  
xticks(1:length(whatDiseases)); yticks(1:length(whatDiseases)); 
xticklabels(whatDiseases); yticklabels(whatDiseases)


subplot(1,2,2); imagesc(r); colormap(colors); caxis([-1 1]); colorbar; 
title('Correlation between coexp'); axis('square'); 
xticks(1:length(whatDiseases)); yticks(1:length(whatDiseases)); 
xticklabels(whatDiseases); yticklabels(whatDiseases)

