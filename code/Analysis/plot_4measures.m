% matching for psychiatric disorders
clear all; close all; 
similarityTypes = {'MAGMAdefault', 'PPI_mapped_th600', 'eQTLbrain', 'AllenMeanCoexpMapped'};
similarityTypes_label = {'SNP position', 'PPI network', 'Brain eQTL', 'AHBA'}; 

whatDiseases_GWAS = {'ADHD','MDD2','SCZ','BIP2','DIABETES'}; 
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

whatMeasures = 'allPsych'; 
% plot correlation between different measures for one representative disorder
%f = figure('color','w', 'Position', [300, 300, 2000, 1000]);
for i=1:numGWAS
    
    %ax{i} = subplot(1,numGWAS, i); hold on;
    [f, Mnames{i}, Mnumbers{i}] = correlate_geneMeasures(whatDiseases_GWAS{i}, whatMeasures, true);
    figureName = sprintf('figures/%s_geneMeasures_%s', whatDiseases_GWAS{i}, whatMeasures);
    print(f,figureName,'-dpng','-r300');
end


% does combinig scores improve matches?









%% don't plot these - will say in words, it's a null result
% f = figure('color','w', 'Position', [300, 300, 1500, 400]); 
% % choose more diverging colours
% barColor = BF_getcmap('set5',5);  
% 
% for i=1:numGWAS
%     
%     ax{i} = subplot(1,numGWAS,i); hold on
%     title(sprintf('%s', whatDiseases_GWAS{i}))
%     b = bar(Ptable.(whatDiseases_GWAS{i}).Pvals);
%     ylabel('-log10(P)')
%     set(gca,'FontSize', 14)
%     b.CData = barColor;
%     b.FaceColor = 'flat';
%     
%     ax{i}.XTick = 1:length(similarityTypes)+1;
%     ax{i}.XTickLabel = [similarityTypes_label, 'Combined'];
%     ax{i}.XTickLabelRotation = 45;
%     
% end
% 
% % rescale axes
% linkaxes([ax{:}],'y');
% 
% print(gcf,'figures/comparePvals','-dpng','-r300');

    
    