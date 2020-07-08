%--------------------------------------------------------------------------- 
% This script imports HMAGMA output files and assigns each gene an entrezID
%--------------------------------------------------------------------------- 

Disorders = {'ADHD', 'MDD2', 'SCZ', 'DIABETES', 'BIP2', 'HF', 'AD'}; 
whatANNOT = {'MAGMAdefault', 'Adult_brain', 'Fetal_brain', 'Neuro', 'Astro', 'eQTLpec'}; 

% load gene ID matching file and select genes that have entrezIDs
entrezID = importGENEIDfile('data/GWASlists/BIOMART_geneIDs.txt'); 
entrezID = entrezID(~isnan(entrezID.entrezgene_id),:); 

entrezIDonly = entrezID(:, {'ensembl_gene_id', 'entrezgene_id', 'external_gene_name'}); 

% there are cases where different gene IDs have the same entrez ID and gene name, selet only usique instances
% of entrezID and gene name combinations and - these are the IDs that will
% be used later. 
entrezIDname = entrezID(:, {'entrezgene_id', 'external_gene_name'}); 
[~, INDunique] = unique(entrezIDname, 'rows');
entrezIDonly = entrezIDonly(INDunique,:); 

for D=1:length(Disorders)
    for A=1:length(whatANNOT)
        
        fileName = sprintf('data/GWASlists/GWASgenes/pgc%s_%s_genes.txt', Disorders{D}, whatANNOT{A});
        mapLIST = importHMAGMAoutfile(fileName); 
        
        fprintf('Processing disorder: %s, annotation %s\n', Disorders{D}, whatANNOT{A});
        
        [~, INDmap, indENT] = intersect(mapLIST.GENE, entrezIDonly.ensembl_gene_id, 'stable'); 
        
        % select only overlapping genes from the table
        % in some cases one gene has multiple entrez IDs, probably related
        % to recent updates, hard to tell which one to choose, so we'll
        % keep gene names as well; 
        mapLIST = mapLIST(INDmap,:); 
        ENTREZID = entrezIDonly.entrezgene_id(indENT);
        GENENAME = entrezIDonly.external_gene_name(indENT);
        
        % add to table
        mapLIST = addvars(mapLIST,ENTREZID,GENENAME,'After','GENE');
        % divided 
        pThr = 0.05/size(mapLIST,1); 
        fprintf('%d genes\n', sum(mapLIST.P<pThr));
        
        % there are also some cases where the same ENTREZis and same gene
        % name corresponds to several GENEIDs, if that's the case, keep one

        writetable(mapLIST, sprintf('data/GWASlists/GWASgenes/pgc%s_%s_genes_entrezID.txt', ...
            Disorders{D}, whatANNOT{A}), 'Delimiter','\t'); 
    end
end


