function measureNames = aggregate_names(diseaseResultsP, similarityTypes, ALLmeasures)

% for every non-nan value, get their measure name
names = cell(length(similarityTypes), length(ALLmeasures)); 
for i=1:length(similarityTypes)
    for j=1:length(ALLmeasures)
        names{i,j} = [similarityTypes{i}, '.', ALLmeasures{j}]; 
    end
end

remIND = isnan(diseaseResultsP); 
N = string(names); 
N(remIND) = 'NaN'; 
INDrem = strcmp(N(:), 'NaN');
measureNames = names(~INDrem); 

end