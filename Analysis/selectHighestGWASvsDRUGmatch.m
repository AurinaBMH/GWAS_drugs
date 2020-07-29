% Loop over different versions of the GWAS-gene mapping and select
% combination that gives the strongest match between GWAS and drugs.
clear all; close all; 

% first run all diseases with itself
whatDiseases_GWAS = {'ADHD', 'MDD2', 'SCZ', 'BIP2', 'DIABETES', 'HF'};
whatNull = 'randomDisease';

for d = 1:length(whatDiseases_GWAS)
    % input requires cell
    whatDisease_GWAS{1} =  whatDiseases_GWAS{d};
    compareGWASvsDRUGmatches(whatDisease_GWAS, whatNull); 
    
    figureName = sprintf('figures/GWAS%s_vs_drug%s_%s', whatDisease_GWAS{1}, whatDisease_GWAS{1}, whatNull);
    print(gcf,figureName,'-dpng','-r300');
    
end

% for PPI select all 1 and 2 neighbor measures
whatDiseases_GWAS = {'BIP2'};
whatDiseases_Drug = 'SCZ';
whatNull = 'randomDisease';

compareGWASvsDRUGmatches(whatDiseases_GWAS, whatNull, whatDiseases_Drug); 
figureName = sprintf('figures/GWAS%s_vs_drug%s_%s', whatDiseases_GWAS{1}, whatDiseases_Drug, whatNull);
print(gcf,figureName,'-dpng','-r300');
    