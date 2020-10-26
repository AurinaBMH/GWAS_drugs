function geneStat = evaluate_allenCoexp(allTargetGenes,GWASgenes)

% Human gene coexpression (AHBA)
%-------------------------------------------------------------------------------
% Look at the distribution of cortical coexpression each gene with the set of
% matches from the list of context genes (e.g., GWAS hits)
%-------------------------------------------------------------------------------
numGenesChar = length(allTargetGenes);
%-------------------------------------------------------------------------------
% Load gene coexpression information from the Allen Human Brain Atlas (processed):
%-------------------------------------------------------------------------------
[geneCoexp,AllenGeneInfo] = LoadCoexpression();

% Find matches:
isGWAS = ismember(AllenGeneInfo.GeneSymbol,GWASgenes);
fprintf(1,'%u/%u context genes could be matched to Allen data\n',sum(isGWAS),length(GWASgenes));

fprintf(1,'%u/%u genes to be characterized successfully matched to Allen data\n',...
                sum(ismember(AllenGeneInfo.GeneSymbol,allTargetGenes)),length(allTargetGenes));

%-------------------------------------------------------------------------------
% Loop over list of genes to characterize:
zval = nan(numGenesChar,1);
log10p = nan(numGenesChar,1);
invp = nan(numGenesChar,1);
mean_r = nan(numGenesChar,1);

for i = 1:numGenesChar
    gene_i = allTargetGenes{i};

    allenIndex = strcmp(AllenGeneInfo.GeneSymbol,gene_i);
    if any(allenIndex)
        % Compute the coexpression values of this gene to the set of context genes:
        % take the absolute values of coexpression
        % all values for GWAS genes
        all_GWAS = abs(geneCoexp(allenIndex,isGWAS));
        % all values for non-GWAS genes
        all_nonGWAS = abs(geneCoexp(allenIndex,~isGWAS));
        
        % evaluate the difference between two distributions - expectation
        % that GWAS genes will have more positive values than non-GWAS genes
        [p,~,stats] = ranksum(all_GWAS, all_nonGWAS, 'tail', 'right'); 
        
        zval(i) = stats.zval; 
        log10p(i) = -log10(p); 
        invp(i) = 1/p; 

        mean_r(i) = nanmean(all_GWAS);

    else
        % This gene could not be matched to AHBA data
        % warning('%s could not be matched to the Allen expression data',gene_i)
        % AllenMeanCoexp(i) = NaN;
    end
end

% save to struct
geneStat.zval = zval; 
geneStat.log10p = log10p; 
geneStat.invp = invp; 
geneStat.mean_r = mean_r; 

end
