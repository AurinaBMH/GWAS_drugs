% matching for psychiatric disorders
clear all; close all; 
similarityTypes = {'MAGMAdefault', 'PPI_mapped_th600', 'eQTLbrain', 'AllenMeanCoexpMapped'}; 
whatDiseases_GWAS = {'ADHD','MDD2','SCZ','BIP2','DIABETES'}; %'BIPandSCZ'

numGWAS = length(whatDiseases_GWAS); 


for s=1:length(similarityTypes)
    
    if contains(similarityTypes{s},'PPI')
        whatProperty = 'percPPIneighbors1';
    else
        if ~contains(similarityTypes{s},'Allen')
            whatProperty = 'P';
        elseif contains(similarityTypes{s},'Allen')
            whatProperty = 'r';
        end
    end
    
    [rhosALL ,pValsALL] = DistinguishingCharBar(similarityTypes{s},whatProperty, 'randomDrugP', 'BF', whatDiseases_GWAS, true, length(similarityTypes));
    figureName = sprintf('figures/BarChart_%s', similarityTypes{s});
    print(gcf,figureName,'-dpng','-r300');
    
    
end

% does combinig scores improve matches?
[Ptable] = compare_optimizedScores(whatDiseases_GWAS, 'reduced', false); 

f = figure('color','w', 'Position', [300, 300, 1500, 400]); 
barColor = [31,120,180;
            33,102,172; 
            233,163,201;
            252,141,89; 
            178,24,43]; 

for i=1:numGWAS
    ax{i} = subplot(1,numGWAS,i); hold on
    b = bar(Ptable.(whatDiseases_GWAS{i}).Pvals);
    
    b.CData = barColor; 
    b.FaceColor = 'flat';
    
end


    
    