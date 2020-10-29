clear all; close all; 

whatDiseases_GWAS = {'ADHD','MDD2','SCZ','BIP2','DIABETES'}; 
similarity = 'AllenMapped'; 
whatMeasure = 'log10p'; 
r = zeros(length(whatDiseases_GWAS)); 
figure('color', 'w');

for i=1:length(whatDiseases_GWAS) 
    load(sprintf('resultsTable_%s_BF_2020_all_drugbank.mat', whatDiseases_GWAS{i}))
    %A1 = 1./10.^-(geneScores.MAGMAdefault.P); 
    A1 = geneScores.(similarity).(whatMeasure); 
    % plot the histogram
    subplot(2,round(length(whatDiseases_GWAS)/2), i); 
    histogram(A1, 20); 
    title(sprintf('%s GWAS', whatDiseases_GWAS{i}))
    xlabel('Gene score')
    
    for j=1:length(whatDiseases_GWAS)
    load(sprintf('resultsTable_%s_BF_2020_all_drugbank.mat', whatDiseases_GWAS{j}))   
        A2 = geneScores.(similarity).(whatMeasure); 
        
        r(i,j) = corr(A1, A2, 'rows', 'complete'); 
        
    end
end

[COL]=cbrewer('div', 'RdBu', 100);
colors = flipud(COL);
figure('color', 'w'); imagesc(r); colormap(colors); 
caxis([-1 1])
xticks(1:length(whatDiseases_GWAS) ); xticklabels(whatDiseases_GWAS); 
yticks(1:length(whatDiseases_GWAS) ); yticklabels(whatDiseases_GWAS); 


