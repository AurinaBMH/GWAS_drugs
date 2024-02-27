function allDrugs = get_allDrugBank_targets(whatTargets)

if nargin < 1
    whatTargets = 'all'; 
end

% this function will get all drugs and their targets from
% pharmacologically_active.csv file and aggregate those to allDrugs format
switch whatTargets
    case 'active'
        targetsBANK = readtable('data/TREATMENTlists/Drug_Bank_database/2024/drugbank_approved_target_polypeptide_ids/pharmacologically_active.csv');
    case 'all'
        targetsBANK = readtable('data/TREATMENTlists/Drug_Bank_database/2024/drugbank_approved_target_polypeptide_ids/all.csv');
end
vocabularyBANK = readtable('data/TREATMENTlists/Drug_Bank_database/2024/drugbank_vocabulary.csv');

% keep only Human genes
IND_human = contains(targetsBANK.Species, 'Humans'); 
targetsBANK = targetsBANK(IND_human,:); 

% find all drugs mentioned in the active target list
activeDrugs = strjoin(targetsBANK.DrugIDs','; ');
activeDrugs = unique(split(activeDrugs, '; '));

% for all drugs find their targets
Name = cell(length(activeDrugs),1); 
Targets = cell(length(activeDrugs),1); 
for i = 1:length(activeDrugs)
    
    % find drug targets
    INDtarg = contains(targetsBANK.DrugIDs, activeDrugs{i});
    targets = targetsBANK.GeneName(INDtarg);
    % find drug name
    INDvoc = find(strcmpi(vocabularyBANK.DrugBankID, activeDrugs{i}), 1);
    
    % remove
    Name{i} = lower(erase(vocabularyBANK.CommonName(INDvoc), {'-', '+', '/','\', ')', '('})); 
    Targets{i} = strjoin(targets,', ');
    
end

allDrugs = table;
allDrugs.Name = Name; 
allDrugs.Target = Targets; 

% remove drugs with non-existent targets
allDrugs(strcmp(string(allDrugs.Target), ""),:) = [];
% select only unique drug names
[~, INDunique] = unique(string(allDrugs.Name));
allDrugs = allDrugs(INDunique,:);


end