function T = give_drugTargets_drugBank(drugName)
% this function finds drug targets for selected drug based on drugBank
% database

% load
vocabulary = readtable('data/TREATMENTlists/Drug_Bank_database/drugbank_vocabulary.csv');
targets = readtable('data/TREATMENTlists/Drug_Bank_database/drugbank_all_target_polypeptide_ids.csv/all.csv');

INDvoc = find(contains(lower(vocabulary.CommonName), lower(drugName)), 1);
if isempty(INDvoc)
    INDvoc = find(contains(lower(vocabulary.Synonyms), lower(drugName)), 1);
end

T = [];

% if found in at least one case
if ~isempty(INDvoc)
    % take drug ID
    drugID = vocabulary.DrugBankID{INDvoc};
    
    % find this drug in another dataset
    INDtarg = contains(targets.DrugIDs, drugID);
    T = unique(targets.GeneName(INDtarg));
end

end