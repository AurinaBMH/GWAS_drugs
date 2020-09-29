function [RNADOMdrugs, geneNames] = give_randomDrug_null(diseaseName, disorderDrugs, allDrugs, whatSelection)

if nargin<4
    whatSelection = 'proportional'; 
end

% This function for a selected disease: 
% 1. finds a number of drugs for that disease
% 2. selects a random set of drugs (same size)
% 3. gives a normalisd score vector for this random set of drugs treating
% it as an additional disease

   
% e.g. for ADHD, it selects 18 random drugs and so on. 
params = SetDefaultParams();
% find the number of drugs selected for that disease
numDrugs = size(disorderDrugs.(diseaseName),1); 

switch whatSelection
    case 'random'
        % now select the same number of drugs at random from the list without replacement
        INDrand = datasample(1:size(allDrugs,1), numDrugs,'Replace',false);
    case 'proportional'
        drugPROB = zeros(size(allDrugs,1),1);
        numDisorders = length(fields(disorderDrugs));
        whatDisorders = fields(disorderDrugs);
        % for each drug, make a probability value of being selected as a
        % proportion of total number of drugs for that disorder
        % this way the probability to select a drug from a particular
        % disorder is equal, so the selected drugs are not over-sampling
        % disorders with long lists of drugs; 
        
        % for each drug, find how many times it's mentioned:
        % 1. if only once, give a score 1/numdisorderDrugs for that disorder
        % 2. if more than once, take a sum of 1/NumdisorderDrugs
        for dr=1:size(allDrugs,1)
            % find whihc disorders a drug is mentioned in
            D = allDrugs.Name{dr};
            idx = zeros(numDisorders,1);
            for p=1:numDisorders
                K = find(ismember(disorderDrugs.(whatDisorders{p}).Name, D));
                if ~isempty(K)
                    idx(p) = K;
                end
            end
            
            % make a weight using 1/numDrugs in disorder
            whatDis = find(idx);
            dpALL = 0;
            for dp = 1:length(whatDis)
                dpSINGLE = 1/size(disorderDrugs.(whatDisorders{whatDis(dp)}),1);
                dpALL = dpALL+dpSINGLE;
            end
            drugPROB(dr) = dpALL;
            
        end
        
        % repeat until all drugs are taken only once - to make it without
        % replacement
        INDrand = datasample(1:size(allDrugs,1), numDrugs,'Replace',false,'Weights',drugPROB); 
    case 'proportionalPsych'
        % find only psychiatry drugs
        allDrugs = vertcat(disorderDrugs.ADHD, disorderDrugs.BIP, disorderDrugs.MDD, disorderDrugs.SCZ);
        [~, ix] = unique(allDrugs.Name, 'stable'); 
        allDrugs = allDrugs(ix,:); 
        
        drugPROB = zeros(size(allDrugs,1),1);
        whatDisorders = {'ADHD'; 'BIP'; 'MDD'; 'SCZ'}; 
        numDisorders = length(whatDisorders);
        
        % for each drug, make a probability value of being selected as a
        % proportion of total number of drugs for that disorder
        % this way the probability to select a drug from a particular
        % disorder is equal, so the selected drugs are not over-sampling
        % disorders with long lists of drugs; 
        
        % for each drug, find how many times it's mentioned:
        % 1. if only once, give a score 1/numdisorderDrugs for that disorder
        % 2. if more than once, take a sum of 1/NumdisorderDrugs
        for dr=1:size(allDrugs,1)
            % find whihc disorders a drug is mentioned in
            D = allDrugs.Name{dr};
            idx = zeros(numDisorders,1);
            for p=1:numDisorders
                K = find(ismember(disorderDrugs.(whatDisorders{p}).Name, D));
                if ~isempty(K)
                    idx(p) = K;
                end
            end
            
            % make a weight using 1/numDrugs in disorder
            whatDis = find(idx);
            dpALL = 0;
            for dp = 1:length(whatDis)
                dpSINGLE = 1/size(disorderDrugs.(whatDisorders{whatDis(dp)}),1);
                dpALL = dpALL+dpSINGLE;
            end
            drugPROB(dr) = dpALL;
            
        end
        
        % repeat until all drugs are taken only once - to make it without
        % replacement
        INDrand = datasample(1:size(allDrugs,1), numDrugs,'Replace',false,'Weights',drugPROB); 
 
end
% make a table for those selected drugs
drugs_rand = allDrugs(INDrand,:); 

normalizeWithinDrugs = true;
% this is a modified version of the original ImportTreatmentLists that adds random set as a separate disease
% need to keep other dseases as well, so the list of targets is complete
indicatorTable = ImportTreatmentLists_random(normalizeWithinDrugs, drugs_rand, params.whatDrugTargets);
geneWeights = indicatorTable.('RANDOM');

% Normalize (non-NaN elements) to unit vector as 2-norm:
whatNorm = 2; 
RNADOMdrugs = normalizeScoreVector(geneWeights, whatNorm); 
geneNames = indicatorTable.Row; 

end

