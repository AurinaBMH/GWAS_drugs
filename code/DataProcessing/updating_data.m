% Process all data to analysis (adding more GWASs, drugs lists and mapping methods)
% 10/09/2020

%1. get targets for drugs; 
dataTable = give_drugTargets('active'); 

%2. make MAGMA file for corresponsing disorders
get_HMAGMAentrezIDs.m

%3. generate results tables for each disorder and each mapping method
GenerateResultsTables

% make new set of nulls for randomDrug
generate_randomDrug_nulls();

% make a modified set of nulls for weighted randomDrug

