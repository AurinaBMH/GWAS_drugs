function T = give_drugTargets_drugBank(vocabulary, targets, drugName)
% this function finds drug targets for selected drug based on drugBank
% database


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