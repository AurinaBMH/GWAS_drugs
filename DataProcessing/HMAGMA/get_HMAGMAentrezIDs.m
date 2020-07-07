%--------------------------------------------------------------------------- 
% This script imports HMAGMA output files and assigns each gene an entrezID
%--------------------------------------------------------------------------- 

Disorders = {'ADHD', 'MDD2', 'SCZ', 'DIABETES', 'BIP2', 'HF', 'AD'}; 
whatANNOT = {'MAGMAdefault', 'Adult_brain', 'Fetal_brain', 'Neuro', 'Astro'}; 

% load gene ID matching file and select genes that have entrezIDs
entrezID = importENTREZIDfile('data/GWASlists/BIOMART_entrezID.txt'); 
entrezID = entrezID(~isnan(entrezID.EntrezID),:); 

for D=1:length(Disorders)
    for A=1:length(whatANNOT)
        
        fileName = sprintf('data/GWASlists/GWASgenes/pgc%s_%s_genes.txt', Disorders{D}, whatANNOT{A});
        mapLIST = importHMAGMAoutfile(fileName); 
        
        [~, INDmap, indENT] = intersect(mapLIST.GENE, entrezID.Gene, 'stable'); 
        
        % select only overlapping genes from the table
        mapLIST = mapLIST(INDmap,:); 
        ENTREZID = entrezID.EntrezID(indENT);
        
        % add to table
        mapLIST = addvars(mapLIST,ENTREZID,'After','GENE');
        
        fprintf('%d genes\n', sum(mapLIST.P<5*10^-6));
        

        writetable(mapLIST, sprintf('data/GWASlists/GWASgenes/pgc%s_%s_genes_entrezID.txt', ...
            Disorders{D}, whatANNOT{A}), 'Delimiter','\t'); 
    end
end


