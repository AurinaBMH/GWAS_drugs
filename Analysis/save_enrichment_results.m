% extract enrichment results
disorders = {'ADHD', 'BIP', 'MDD', 'SCZ'}; 
measureType = 'PPI'; 

for d=1:length(disorders)
    
GOtableCON = makeEnrichmentResultTable(disorders{d}, measureType); 

end
