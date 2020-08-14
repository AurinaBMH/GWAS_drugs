function plot_drugScores()

whatDiseases_Treatment = {'ADHD','BIP','SCZ','MDD','pulmonary','cardiology','gastro','diabetes'};
geneWeightsNormDrug = cell(length(whatDiseases_Treatment),1); 

for i=1:length(whatDiseases_Treatment)
    
    [~,geneWeightsNormDrug{i}] = GiveMeNormalizedScoreVector(whatDiseases_Treatment{i},'Drug');
    
end

A = cat(2, geneWeightsNormDrug{:}); 
figure('color','w'); imagesc(A); 
ylabel('Genes')
xticks(1:length(whatDiseases_Treatment))
xticklabels(whatDiseases_Treatment); 
colormap(cbrewer('seq', 'Blues', 64));
caxis([0 0.1])

end
