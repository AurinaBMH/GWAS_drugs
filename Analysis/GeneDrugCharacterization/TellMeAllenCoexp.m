function AllenMeanCoexp = TellMeAllenCoexp(genesChar,genesContext)
% Human gene coexpression (AHBA)
%-------------------------------------------------------------------------------
% Look at the distribution of cortical coexpression each gene with the set of
% matches from the list of context genes (e.g., GWAS hits)
%-------------------------------------------------------------------------------
numGenesChar = length(genesChar);

%-------------------------------------------------------------------------------
% Load gene coexpression information from the Allen Human Brain Atlas (processed):
%-------------------------------------------------------------------------------
[geneCoexp,AllenGeneInfo] = LoadCoexpression();

% Find matches:
isContext = ismember(AllenGeneInfo.GeneSymbol,genesContext);
fprintf(1,'%u/%u context genes could be matched to Allen data\n',sum(isContext),length(genesContext));

fprintf(1,'%u/%u genes to be characterized successfully matched to Allen data\n',...
                sum(ismember(AllenGeneInfo.GeneSymbol,genesChar)),length(genesChar));

%-------------------------------------------------------------------------------
% Loop over list of genes to characterize:
AllenMeanCoexp = nan(numGenesChar,1);

for i = 1:numGenesChar
    gene_i = genesChar{i};

    allenIndex = strcmp(AllenGeneInfo.GeneSymbol,gene_i);
    if any(allenIndex)
        % Compute the coexpression values of this gene to the set of context genes:
        coExpContext = geneCoexp(allenIndex,isContext);
        AllenMeanCoexp(i) = nanmean(coExpContext);
    else
        % This gene could not be matched to AHBA data
        % warning('%s could not be matched to the Allen expression data',gene_i)
        % AllenMeanCoexp(i) = NaN;
    end
end

end
