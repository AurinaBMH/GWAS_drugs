% extract enrichment results
function GOtable = save_enrichment_results(measureType)
if nargin <1
    measureType = 'PPI';
end

disorders = {'ADHD', 'BIP', 'MDD', 'SCZ'};

for d=1:length(disorders)
    
    GOtable.(disorders{d}) = makeEnrichmentResultTable(disorders{d}, measureType);
    
end
end
