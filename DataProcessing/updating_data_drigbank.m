% reprocessing data using DrugBank data: 
dataTable = give_drugTargets('active', 'drugbank'); 
GenerateResultsTables

% generate nulls
generate_randomDrug_nulls('drugbank')

% generate results

