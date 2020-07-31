function data = give_drugTargets()
% This function gives drug targets for selected lists of disorders based on
% DrugBank and Drug repurposing hub
disorders = {'ADHD', 'BIP', 'cardiology', 'diabetes', 'gastro', 'MDD', 'pulmonary', 'SCZ'};
% import DrugBank database

% import Drug repurposing hub
drugREP = importDrug_rep_hub('data/TREATMENTlists/Drug_repurposing_hub_database/Drug_repurposing_hub_database_20200730.txt'); 
data = struct; 

for d = 1:length(disorders)
    % treatment_list_%s.txt contain treatments used for each disorder
    
    listDrugs = readcell(sprintf('data/TREATMENTlists/Drug_list_disorders/treatment_list_%s.txt', disorders{d}));
    
    for i = 1:length(listDrugs)
        [~, INDrep] = intersect(drugREP.pert_iname, listDrugs{i}); 
        [~, INDbank] = intersect(drugBANK.XXX, listDrugs{i}); 
        % check in one database
        if ~isempty(INDrep)
        targetREP = drugREP.target(INDrep);     
        targetREPlist = strsplit(targetREP, '|');   
        end
        % check in another database
        if ~isempty(INDbank)
        targetBANK = drugBANK.target(INDbank);     
        targetBANKlist = strsplit(targetBANK, '|');   
            
        end
        
        targetCOMB = unique(horzcat(targetREPlist, targetBANKlist)); 
        data.(disorders{d}).(listDrugs{i}) = targetCOMB'; 
        
    end
    
    
end

end