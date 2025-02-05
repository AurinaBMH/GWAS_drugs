% Generating results for psychiatric disorders:
function Pmatrix = plot_sensitivity_figures(whatDisorder)

if nargin < 1
    whatDisorder = 'BIP'; 
end
%-------------------------------------------------------
% Set the options
%-------------------------------------------------------
params = SetDefaultParams();
similarityTypes = {'MAGMAdefault', 'PPI_mapped_th600', 'eQTLbrain', 'AlleneQTLbrain'};
similarityTypes_label = {'SNP position', 'PPI network', 'Brain eQTL', 'AHBA'}; 

switch whatDisorder
    case 'BIP'
        treatment_classes = params.whatDiseases_Treatment_classes; 
    case 'MDD'
        treatment_classes = params.whatDiseases_Treatment_classes_MDD;
    case 'ADHD'
        treatment_classes = params.whatDiseases_Treatment_classes_ADHD;
    case 'SCZ'
        treatment_classes = params.whatDiseases_Treatment_classes_SCZ;
end

ix_dis = contains(params.whatGWAS_2024, whatDisorder);
whatDisease_GWAS = params.whatGWAS_2024(ix_dis);

numDrugs = length(treatment_classes);
whatMeasures = 'treatment_class';
whatNull = sprintf('randomDrugR_%s_drugbank_treatment_class', params.whatTargets);

%-------------------------------------------------------
% Figure S1: Calculate matches between BIP3 GWAS and each treatment class
%-------------------------------------------------------

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
    
    [rhosALL, pValsALL_randDrug, whatDiseases_Treatment_randDrug, nullScoresALL{s}] = DistinguishingCharBar_sensitivity(similarityTypes{s}, whatProperty, whatNull, 'BF',...
        whatDisease_GWAS, false, numDrugs, whatMeasures, whatDisorder);

        
    % find corresponsing match
    [T, INDr] = intersect(whatDiseases_Treatment_randDrug, treatment_classes, 'stable'); 
    % select disorder to itself - diagonal
    Pmatrix(s,:) = pValsALL_randDrug(INDr); 
    RHOmatrix(s,:) = rhosALL(INDr); 

end

%-------------------------------------------------------
% Figure 2: Calculate matches for individual disorders
%-------------------------------------------------------
f = plot_measureOverview_sensitivity(Pmatrix, T, similarityTypes_label, whatDisorder); 
figureName = sprintf('figures_2024/BarP_%s_treatment_class_%s_%s', whatDisorder, whatMeasures, whatNull);
print(gcf,figureName,'-dpng','-r300');


end


    
    