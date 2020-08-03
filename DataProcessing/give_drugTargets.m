function data = give_drugTargets()
% This function gives drug targets for selected lists of disorders based on
% DrugBank and Drug repurposing hub
disorders = {'ADHD', 'BIP', 'cardiology', 'diabetes', 'gastro', 'MDD', 'pulmonary', 'SCZ'};
% import Drug repurposing hub
drugREP = importDrug_rep_hub('data/TREATMENTlists/Drug_repurposing_hub_database/Drug_repurposing_hub_database_20200730.txt');

% load drugBank database files
vocabularyBANK = readtable('data/TREATMENTlists/Drug_Bank_database/drugbank_vocabulary.csv');
%targetsBANK = readtable('data/TREATMENTlists/Drug_Bank_database/drugbank_all_target_polypeptide_ids.csv/pharmacologically_active.csv');
targetsBANK = readtable('data/TREATMENTlists/Drug_Bank_database/drugbank_all_target_polypeptide_ids.csv/all.csv');
data = struct;

for d = 1:length(disorders)
    
    % treatment_list_%s.txt contain treatments used for each disorder
    % files are made manually based on drug lists that Jannette selected
    % these are just the first columns from, e.g.
    % Treatment-list-ADHD-4thMay2018.xlsx + a couple of additional drugs
    % from Ken's email on the 29/07/2020: 
    % ADHD: +mydayis
    % SCZ:  +lumateperone
    % MDD:  +esketamine and brexanolone
    
    listDrugs = readcell(sprintf('data/TREATMENTlists/Drug_list_disorders/treatment_list_%s.txt', disorders{d}));
    
    for i = 1:length(listDrugs)
        
        % get targets from hub
        [~, INDrep] = intersect(drugREP.pert_iname, listDrugs{i});
        % check in one database
        if ~isempty(INDrep)
            
            targetREP = drugREP.target(INDrep);
            targetREPlist = strsplit(targetREP, '|')';
        else
            
            targetREPlist = string(); 
        end
        % get targets from Drug Bank
        targetBANKlist = string(give_drugTargets_drugBank(vocabularyBANK, targetsBANK, listDrugs{i})); 
        
        % combine both and select unique genes
        targetCOMB = vertcat(targetREPlist, targetBANKlist);
        % in case there are differences in cases, make all upper and only the select unique
        [~, ia] = unique(upper(targetCOMB));
        targetCOMB = targetCOMB(ia);
        
        % remove some characters to make it a viable structure field
        drugName = erase(listDrugs{i}, {'-', '+', '/','\', ')', '('}); 
        data.(disorders{d}).(drugName) = targetCOMB'; 
        
    end
    
    
end

end