
whatDiseases_GWAS = {'ADHD','MDD2','SCZ','BIP2','DIABETES'};
r = zeros(length(whatDiseases_GWAS)); 

for i=1:length(whatDiseases_GWAS) 
    load(sprintf('resultsTable_%s_BF_2020_all_drugbank.mat', whatDiseases_GWAS{i}))
    A1 = geneScores.PPI_mapped_th600.percPPIneighbors1;
    for j=1:length(whatDiseases_GWAS)
    load(sprintf('resultsTable_%s_BF_2020_all_drugbank.mat', whatDiseases_GWAS{j}))   
        A2 = geneScores.PPI_mapped_th600.percPPIneighbors1;
        
        r(i,j) = corr(A1, A2, 'rows', 'complete'); 
        
    end
end
[COL]=cbrewer('div', 'RdBu', 100);
colors = flipud(COL);
figure; imagesc(r); colormap(colors); 
caxis([-1 1])
xticks(1:5); xticklabels(whatDiseases_GWAS); 
yticks(1:5); yticklabels(whatDiseases_GWAS); 