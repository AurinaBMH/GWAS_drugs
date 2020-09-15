% Loop over different versions of the GWAS-gene mapping and select
% combination that gives the strongest match between GWAS and drugs.
clear all; close all;

params = SetDefaultParams(); 
% first run all diseases with itself
whatDiseases_GWAS = params.whatGWAS; 
whatNull = 'randomDrugP';

for dg = 1:length(whatDiseases_GWAS)
    
    whatDisease_GWAS{1} =  whatDiseases_GWAS{dg};
    
    compareGWASvsDRUGmatches(whatDisease_GWAS, whatNull);
    
    % this is for one
    figureName = sprintf('figures/GWAS%s_vs_drug_%s_%stargets', whatDisease_GWAS{1}, whatNull, params.whatDrugTargets);
    print(gcf,figureName,'-dpng','-r300');
    
end


% run pairs of diseases

whatDiseases_GWAS = {'MDD2', 'SCZ', 'BIP2'}; 
whatDiseases_Treatment = {'BIP','SCZ','MDD'};
whatNull = 'randomDrugP';

diseaseResultsR = cell(length(whatDiseases_GWAS), length(whatDiseases_Treatment));
diseaseResultsP = cell(length(whatDiseases_GWAS), length(whatDiseases_Treatment));

for dg = 1:length(whatDiseases_GWAS)
  for dt = 1:length(whatDiseases_Treatment)
    % input requires cell
    whatDisease_GWAS{1} =  whatDiseases_GWAS{dg};
    Dname = whatDiseases_Treatment{dt};
    [diseaseResultsR{dg, dt}, diseaseResultsP{dg, dt}] = compareGWASvsDRUGmatches(whatDisease_GWAS, whatNull, Dname);
    
    % this is for one
    figureName = sprintf('figures/GWAS%s_vs_drug%s_%s', whatDisease_GWAS{1}, Dname, whatNull);
    print(gcf,figureName,'-dpng','-r300');
    
  end
end


% rows are GWAS lists, columns are drugs
fileName = sprintf('DataOutput/GWASvsDRUGS_%s.mat', whatNull); 
save(fileName, 'diseaseResultsR', 'diseaseResultsP', 'whatDiseases_GWAS', 'whatDiseases_Treatment')



% for PPI select all 1 and 2 neighbor measures
% whatDiseases_GWAS = {'BIP2'};
% whatDiseases_Drug = 'MDD2';
% whatNull = 'randomDrug';
% 
% compareGWASvsDRUGmatches(whatDiseases_GWAS, whatNull, whatDiseases_Drug);
% figureName = sprintf('figures/GWAS%s_vs_drug%s_%s', whatDiseases_GWAS{1}, whatDiseases_Drug, whatNull);
% print(gcf,figureName,'-dpng','-r300');
