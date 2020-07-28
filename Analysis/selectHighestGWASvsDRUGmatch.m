% Loop over different versions of the GWAS-gene mapping and select
% combination that gives the strongest match between GWAS and drugs.

load('resultsTable_ADHD_FDR.mat', 'geneScores')
% select all available similarity types
similarityTypes = setdiff(fieldnames(geneScores), {'gene', 'params'});
% for PPI select 1 neighbor and mean and median measures;
PPImeasures = contains(fieldnames(geneScores.PPI_eQTLbrain_th0), '1') | ...
    contains(fieldnames(geneScores.PPI_eQTLbrain_th0), 'mean') | ...
    contains(fieldnames(geneScores.PPI_eQTLbrain_th0), 'median');
PPInum = find(PPImeasures);
whatThreshold = 'BF';
whatNull = 'randomDisease';
whatDiseases_GWAS = {'BIP2'};

diseaseResultsR = nan(length(similarityTypes), length(PPInum));
diseaseResultsP = nan(length(similarityTypes), length(PPInum));

Dname = whatDiseases_GWAS{1};
Dname = Dname(isstrprop(Dname,'alpha'));

for s=1:length(similarityTypes)
    
    if contains(similarityTypes{s},'PPI')
        for p=1:length(PPInum)
            whatProperties = fieldnames(geneScores.(similarityTypes{s}));
            whatProperty = whatProperties{PPInum(p)};
            [rhos ,pVals, whatDiseases_Treatment] = DistinguishingCharBar(similarityTypes{s},whatProperty, whatNull, whatThreshold, whatDiseases_GWAS, false);
            % select rho and p values for a selected disorder
            
            takeVal = contains(whatDiseases_Treatment, Dname, 'IgnoreCase',true);
            
            diseaseResultsR(s,p) = rhos(takeVal);
            diseaseResultsP(s,p) = pVals(takeVal);
            
        end
    elseif contains(similarityTypes{s},'Allen')
        whatProperty = 'r';
        [rhos ,pVals, whatDiseases_Treatment] = DistinguishingCharBar(similarityTypes{s},whatProperty, whatNull, whatThreshold, whatDiseases_GWAS, false);
        takeVal = contains(whatDiseases_Treatment, Dname, 'IgnoreCase',true);
        
        diseaseResultsR(s,1) = rhos(takeVal);
        diseaseResultsP(s,1) = pVals(takeVal);
        
        
    else
        whatProperty = 'P';
        [rhos ,pVals, whatDiseases_Treatment] = DistinguishingCharBar(similarityTypes{s},whatProperty, whatNull, whatThreshold, whatDiseases_GWAS, false);
        takeVal = contains(whatDiseases_Treatment, Dname, 'IgnoreCase',true);
        
        diseaseResultsR(s,1) = rhos(takeVal);
        diseaseResultsP(s,1) = pVals(takeVal);
        
    end
end

figure; imagesc(diseaseResultsR); 
hold on; 
plot_matrixValues(diseaseResultsP)



