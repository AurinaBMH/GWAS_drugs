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
whatMeasures = 'allPsych';
whatNull = sprintf('randomDrugR_%s_drugbank_treatment_class', params.whatTargets); 
numGWAS = length(whatDiseases_classes); 
V = nan(length(similarityTypes), numGWAS); 
whatDiseases_GWAS_name = cell(length(whatDiseases_classes), 1); 

% give names without numbers
for i=1:length(whatDiseases_classes) 
    name = whatDiseases_classes{i};
    whatDiseases_GWAS_name{i} = name(isstrprop(name,'alpha')); 
end

%-------------------------------------------------------
% Figure S1: Calculate pairwise matches 
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
    
    [rhosALL, pValsALL_randDrug, whatDiseases_Treatment_randDrug] = DistinguishingCharBar_sensitivity(similarityTypes{s},whatProperty, whatNull, 'BF', whatDiseases_classes, true, numDrugs, whatMeasures);

        
    % find corresponsing match
    [T, INDr, INDc] = intersect(whatDiseases_Treatment_randDrug, whatDiseases_GWAS_name, 'stable'); 
    % select disorder to itself - diagonal
    Pmatrix(s,:) = diag(pValsALL_randDrug(INDr, INDc)); 

    figureName = sprintf('figures_2024/BarChart_psych_%s_%s_%s', similarityTypes{s}, whatMeasures, whatNull);
    print(gcf,figureName,'-dpng','-r300');


%-------------------------------------------------------
% Figure 2: Calculate matches for individual disorders
%-------------------------------------------------------
f = plot_measureOverview(Pmatrix, T, similarityTypes_label); 
figureName = sprintf('figures_2024/BarP_withinDisorder_%s_%s', whatMeasures, whatNull);
print(gcf,figureName,'-dpng','-r300');


%-------------------------------------------------------
% Figure 3: Correspondence across different data processing methods
%-------------------------------------------------------
DOrecalc = true; % select true if any data was updated since the last run
pPlot_all = plot_compareMeasures(whatDiseases_classes, whatMeasures, DOrecalc); 


end


    
    