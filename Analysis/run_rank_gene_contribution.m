clear all; close all; 
%similarityTypes = {'MAGMAdefault', 'PPI_mapped_th600', 'eQTLbrain', 'AllenMeanCoexpMapped'};
% rank genes for relevant pairs of GWAS vs drug on SNP position measure
similarityType = 'MAGMAdefault';
% these lists are pairs of GWASvsDrugs for ranking;
whatGWAS = {'MDD2', 'SCZ', 'DIABETES'};
whatDrugs = {'BIP', 'BIP', 'DIABETES'};
for i=1:length(whatDrugs)
    T_SNPposition.(whatGWAS).(whatDrugs) = rank_gene_contribution(whatGWAS{i}, whatDrugs{i}, similarityType);
end


% rank genes for relevant pairs of GWAS vs drug on PPI network
similarityType = 'PPI_mapped_th600';
% these lists are pairs of GWASvsDrugs for ranking;
whatGWAS = {'SCZ', 'BIP2', 'BIP2', 'DIABETES'};
whatDrugs = {'BIP', 'BIP', 'SCZ', 'DIABETES'}; 

for i=1:length(whatDrugs)
    T_PPInetwork.(whatGWAS).(whatDrugs) = rank_gene_contribution(whatGWAS{i}, whatDrugs{i}, similarityType);
end



