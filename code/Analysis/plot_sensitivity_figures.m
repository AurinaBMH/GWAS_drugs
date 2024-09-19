% Generating results for psychiatric disorders:
function Pmatrix = plot_sensitivity_figures()
%-------------------------------------------------------
% Set the options
%-------------------------------------------------------
params = SetDefaultParams();
similarityTypes = {'MAGMAdefault', 'PPI_mapped_th600', 'eQTLbrain', 'AlleneQTLbrain'};
similarityTypes_label = {'SNP position', 'PPI network', 'Brain eQTL', 'AHBA'}; 
whatDiseases_classes = params.whatDiseases_Treatment_classes;
numDrugs = length(whatDiseases_classes); 
whatMeasures = 'BIP_treatment_class';
whatNull = sprintf('randomDrugR_%s_drugbank_treatment_class', params.whatTargets); 

%-------------------------------------------------------
% Figure S1: Calculate matches between BIP3 GWAS and each treatment class
%-------------------------------------------------------

nullScores = cell(length(similarityTypes),1); 

for s=1:length(similarityTypes)
    
    if contains(similarityTypes{s},'PPI')
        whatProperty = 'percPPIneighbors1';
    else
        if ~contains(similarityTypes{s},'Allen')
            whatProperty = 'P';
        elseif contains(similarityTypes{s},'Allen')
            whatProperty = 'zval';
        end
    end
    
    [rhosALL, pValsALL_randDrug, whatDiseases_Treatment_randDrug] = DistinguishingCharBar_sensitivity(similarityTypes{s}, whatProperty, whatNull, 'BF', {'BIP3'}, true, numDrugs, whatMeasures);

        
    % find corresponsing match
    [T, INDr] = intersect(whatDiseases_Treatment_randDrug, whatDiseases_classes, 'stable'); 
    % select disorder to itself - diagonal
    Pmatrix(s,:) = pValsALL_randDrug(INDr); 

    figureName = sprintf('figures_2024/BarChart_BIPsensitivity_%s_%s_%s', similarityTypes{s}, whatMeasures, whatNull);
    print(gcf,figureName,'-dpng','-r300');
end

%-------------------------------------------------------
% Figure 2: Calculate matches for individual disorders
%-------------------------------------------------------
f = plot_measureOverview_sensitivity(Pmatrix, T, similarityTypes_label); 
figureName = sprintf('figures_2024/BarP_withinDisorder_%s_%s', whatMeasures, whatNull);
print(gcf,figureName,'-dpng','-r300');
% 
% 
% %-------------------------------------------------------
% % Figure 3: Correspondence across different data processing methods
% %-------------------------------------------------------
% DOrecalc = true; % select true if any data was updated since the last run
% pPlot_all = plot_compareMeasures(whatDiseases_classes, whatMeasures, DOrecalc); 

end


    
    