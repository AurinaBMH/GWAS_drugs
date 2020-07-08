%--------------------------------------------------------------------------- 
% This script imports HMAGMA output files and assigns each gene an entrezID
%--------------------------------------------------------------------------- 

Disorders = {'ADHD', 'MDD2', 'SCZ', 'DIABETES', 'BIP2', 'HF', 'AD'}; 
whatANNOT = {'MAGMAdefault', 'Adult_brain', 'Fetal_brain', 'Neuro', 'Astro', 'eQTLpec'}; 

% load gene ID matching file and select genes that have entrezIDs
entrezID = importGENEIDfile('data/GWASlists/BIOMART_geneIDs.txt'); 
entrezIDonly = entrezID(:, {'ensembl_gene_id', 'entrezgene_id'}); 

entrezIDonly = entrezIDonly(~isnan(entrezIDonly.entrezgene_id),:); 
entrezIDonly = unique(entrezIDonly, 'rows');

for D=1:length(Disorders)
    for A=1:length(whatANNOT)
        
        fileName = sprintf('data/GWASlists/GWASgenes/pgc%s_%s_genes.txt', Disorders{D}, whatANNOT{A});
        mapLIST = importHMAGMAoutfile(fileName); 
        
        [~, INDmap, indENT] = intersect(mapLIST.GENE, entrezIDonly.ensembl_gene_id, 'stable'); 
        
        % select only overlapping genes from the table
        mapLIST = mapLIST(INDmap,:); 
        ENTREZID = entrezIDonly.entrezgene_id(indENT);
        
        % add to table
        mapLIST = addvars(mapLIST,ENTREZID,'After','GENE');
        % divided 
        pThr = 0.05/size(mapLIST,1); 
        fprintf('%d genes\n', sum(mapLIST.P<pThr));
        

        writetable(mapLIST, sprintf('data/GWASlists/GWASgenes/pgc%s_%s_genes_entrezID.txt', ...
            Disorders{D}, whatANNOT{A}), 'Delimiter','\t'); 
    end
end


