function measureNames = aggregate_names(diseaseResultsP, similarityTypes, PPImeasures_names)

% for every non-nan value, get their measure name


k=1;
for s=1:length(similarityTypes)
    if contains(similarityTypes{s},'PPI')
        for p=1:length(PPImeasures_names)
            
            whatProperty = PPImeasures_names{p};
            similarityTypes{s} = strrep(similarityTypes{s},'_',' '); 
            measureNames{k} = [similarityTypes{s}, '.' whatProperty]; 
            k=k+1;
        end
    else
        if ~contains(similarityTypes{s},'Allen')
            whatProperty = 'P';
        elseif contains(similarityTypes{s},'Allen')
            whatProperty = 'r';
        end
        measureNames{k} = [similarityTypes{s}, '.', whatProperty]; 
        k=k+1;
    end
end