function RNADOMdrugs = give_randomDrug_null(diseaseName, disorderDrugs, allDrugs)
% This function for a selected disease: 
% 1. finds a number of drugs for that disease
% 2. selects a random set of drugs (same size)
% 3. gives a normalisd score vector for this random set of drugs treating
% it as an additional disease

% e.g. for ADHD, it selects 18 random drugs and so on. 

% find the number of drugs selected for that disease
numDrugs = size(disorderDrugs.(diseaseName),1); 

% now select the same number of drugs at random from the list without replacement
INDrand = randsample(1:size(allDrugs,1),numDrugs, false); 
% make a table for those selected drugs
drugs_rand = allDrugs(INDrand,:); 

normalizeWithinDrugs = true;
% this is a modified version of the original ImportTreatmentLists that adds random set as a separate disease
% need to keep other idseases as well, so the list of targets is complete
indicatorTable = ImportTreatmentLists_random(normalizeWithinDrugs, drugs_rand);
geneWeights = indicatorTable.('RANDOM');

% Normalize (non-NaN elements) to unit vector as 2-norm:
whatNorm = 2; 
RNADOMdrugs = normalizeScoreVector(geneWeights, whatNorm); 

end

