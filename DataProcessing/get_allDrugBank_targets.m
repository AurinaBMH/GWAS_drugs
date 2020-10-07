function allDrugs = get_allDrugBank_targets()
% this function will get all drugs and their targets from
% pharmacologically_active.csv file and aggregate those to allDrugs format
targetsBANK = readtable('data/TREATMENTlists/Drug_Bank_database/drugbank_all_target_polypeptide_ids.csv/pharmacologically_active.csv');
vocabularyBANK = readtable('data/TREATMENTlists/Drug_Bank_database/drugbank_vocabulary.csv');

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
    
    Name{i} = vocabularyBANK.CommonName(INDvoc);
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