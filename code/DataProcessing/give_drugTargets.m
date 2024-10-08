function dataTable = give_drugTargets(whatTargets, whatDatabase)
if nargin < 1
    whatTargets = 'all'; 
end
if nargin < 2
    whatDatabase = 'drugbank'; 
end

% currently the total list of genes is based on all pharmacologically
% active targets, so including drug repurposing hub information is not
% relevant - these genes will not be used. 

% This function gives drug targets for selected lists of disorders based on
% DrugBank and Drug repurposing hub
params = SetDefaultParams();
disorders = params.whatDiseases_Treatment_ALL;

% import Drug repurposing hub
if strcmp(whatDatabase, 'both')
    drugREP = importDrug_rep_hub('data/TREATMENTlists/Drug_repurposing_hub_database/Drug_repurposing_hub_database_20200730.txt');
    IND_Launched = contains(string(drugREP.clinical_phase), 'Launched');
    drugREP = drugREP(IND_Launched,:);
end
% load drugBank database files
vocabularyBANK = readtable('data/TREATMENTlists/Drug_Bank_database/2024/drugbank_vocabulary.csv');
%
switch whatTargets
    case 'active'
        targetsBANK = readtable('data/TREATMENTlists/Drug_Bank_database/2024/drugbank_approved_target_polypeptide_ids/pharmacologically_active.csv');
    case 'all'
        targetsBANK = readtable('data/TREATMENTlists/Drug_Bank_database/2024/drugbank_approved_target_polypeptide_ids/all.csv');
end

IND_human = contains(targetsBANK.Species, 'Humans'); 
targetsBANK = targetsBANK(IND_human,:); 

dataTable = struct;

for d = 1:length(disorders)
    
    % treatment_list_%s.txt contain treatments used for each disorder
    % files are made manually based on drug lists that Ken selected
    % these are just the first columns from, e.g.
    % Drug lists_KP_September2020.xlsx 
    
    listDrugs = readcell(sprintf('data/TREATMENTlists/Drug_list_disorders/treatment_list_%s.txt', disorders{d}));
    drugName = cell(length(listDrugs),1); 
    targetCOMB = cell(length(listDrugs),1);
    k=1; 
    for i = 1:length(listDrugs)
    
        if strcmp(whatDatabase, 'both')
            % get targets from hub
            [~, INDrep] = intersect(drugREP.pert_iname, listDrugs{i});
            % check in one database
            if ~isempty(INDrep)
                
                targetREP = drugREP.target(INDrep);
                targetREPlist = strsplit(targetREP, '|')';
            else
                
                targetREPlist = string();
            end
        else
            % assign an empty cell, Drug reporposing hub in not relevant;;
            targetREPlist = [];
        end
        % get targets from Drug Bank
        targetBANKlist = string(give_drugTargets_drugBank(vocabularyBANK, targetsBANK, listDrugs{i})); 
        
        % combine both and select unique genes
        targetCOMB_s = vertcat(targetREPlist, targetBANKlist);
                
        % if there are no targets for the drug, skip
        if ~all(strcmp(targetCOMB_s, ""))
        targetCOMB_s(strcmp(targetCOMB_s,"")) = [];
        
        % in case there are differences in cases, make all upper and only the select unique
        [~, ia] = unique(upper(targetCOMB_s));
        targetCOMB_s = targetCOMB_s(ia);
        
        % combine all into one string to be compatible with Ben's old script
        targetCOMB_s = strjoin(targetCOMB_s,',');

        
        % remove some characters to make it a viable structure field
        drugName{k} = erase(listDrugs{i}, {'-', '+', '/','\', ')', '('}); 
        targetCOMB{k} = targetCOMB_s; 
        k=k+1; 
        end

        
    end
    % there are some drugs that have no targets, remove those
    R = find(cellfun(@isempty, drugName)); 
    drugName(R) = [];
    targetCOMB(R) = []; 
    
    % put it in the same format as Ben did before: structure of tables
    dataTable.(disorders{d}) = table(drugName, targetCOMB,'VariableNames',{'Name','Target'});
    fprintf(1,'%s has %u drugs\n',disorders{d},length(drugName));
end

% active is default - will not have a not in the name, all will have a flag
% in the file name

fileName = sprintf('DataOutput_2024/drugTargets_%s_%s_%s.mat', params.whatDrugTargets, whatTargets, whatDatabase);
        
save(fileName, 'dataTable'); 

end