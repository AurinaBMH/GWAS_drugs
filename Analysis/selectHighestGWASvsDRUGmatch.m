% Loop over different versions of the GWAS-gene mapping and select
% combination that gives the strongest match between GWAS and drugs.
clear all; close all;

% first run all diseases with itself
whatDiseases_GWAS = {'ADHD', 'MDD2', 'SCZ', 'BIP2', 'DIABETES', 'HF'};
whatDiseases_Treatment = {'ADHD','BIP','SCZ','MDD','pulmonary','cardiology','gastro','diabetes'};
whatNull = 'randomDrug';

diseaseResultsR = cell(length(whatDiseases_GWAS), length(whatDiseases_Treatment));
diseaseResultsP = cell(length(whatDiseases_GWAS), length(whatDiseases_Treatment));

for d = 1:length(whatDiseases_GWAS)
  for dt = 1:length(whatDiseases_Treatment)
    % input requires cell
    whatDisease_GWAS{1} =  whatDiseases_GWAS{d};
    Dname = whatDiseases_Treatment{dt};
    [diseaseResultsR{d, dt}, diseaseResultsP{d, dt}] = compareGWASvsDRUGmatches(whatDisease_GWAS, whatNull, Dname);

    figureName = sprintf('figures/GWAS%s_vs_drug%s_%s', whatDisease_GWAS{1}, whatDisease_GWAS{1}, whatNull);
    print(gcf,figureName,'-dpng','-r300');
  end
end

% rows are GWAS lists, columns are drugs
fileName = sprintf('DataOutput/GWASvsDRUGS_%s'. whatNull)
save(fileName, 'diseaseResultsR', 'diseaseResultsP', 'whatDiseases_GWAS', 'whatDiseases_Treatment')

% for PPI select all 1 and 2 neighbor measures
whatDiseases_GWAS = {'BIP2'};
whatDiseases_Drug = 'SCZ';
whatNull = 'randomDisease';

compareGWASvsDRUGmatches(whatDiseases_GWAS, whatNull, whatDiseases_Drug);
figureName = sprintf('figures/GWAS%s_vs_drug%s_%s', whatDiseases_GWAS{1}, whatDiseases_Drug, whatNull);
print(gcf,figureName,'-dpng','-r300');
