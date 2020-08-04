function give_randomDrug_null(diseaseName, disorderDrugs, allDrugs)
% for each disease get a random set of drugs that is the same size as 
% real list of drugs, e.g. for SCZ select 45 drugs; 

% find the number of drugs selected for that disease
numDrugs = size(disorderDrugs.(diseaseName),1); 
% now select the same number of drugs at random from the list without replacement
INDrand = randsample(1:size(allDrugs,1),numDrugs, false); 
% make a table for those selected drugs
drugs_rand = allDrugs(INDrand,:); 

end

