% calculate the number of targets per drug
allDrugs = get_allDrugBank_targets('all'); 
load('drugTargets_2020_all_drugbank.mat')

numDis = length(fields(dataTable)); 
dis = fields(dataTable); 

NrT_dis = cell(numDis,1); 
figure('color', 'w'); 
for i=1:numDis
    Drug = dataTable.(dis{i}); 
    NrT = zeros(size(Drug,1),1); 
    for d=1:size(Drug,1)
        T = Drug.Target{d}; 
        NrT(d) = length(strsplit(T, ',')); 
    end
    subplot(2,5,i); histogram(NrT);
    title(sprintf('%s', dis{i})); box off; 
    xlabel('Number of targets per drug')
    ylabel('Number of drugs')
    xlim([0 50])
    NrT_dis{i} = NrT; 
end

% is there a higher proportion of all drugs that have more targets


