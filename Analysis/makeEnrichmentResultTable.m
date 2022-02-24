function GOtableCON = makeEnrichmentResultTable(disorder, measuretype)
% Compile ermineJ results and keep only signifficant values after FDR correction

% this is output file from ermineJ
if strcmp(measuretype, 'PPI')
    analysistype = 'ORA';
    onwhat = {'drugs', 'gwas'};
elseif strcmp(measuretype, 'POS')
    analysistype = 'GSR';
    onwhat = {'gwas'};
end

for n=1:length(onwhat)
    
    type = sprintf('%s_%s_run_on_%s_%s.erminej', disorder, analysistype, onwhat{n}, measuretype);
    fileINname = sprintf('%s.txt', type);
    
    ermineJResultsBPcon = ReadInErmineJ(sprintf('enrichment/output/%s', fileINname));
    CONall = vertcat(ermineJResultsBPcon);
    CONall.corr_pval = mafdr(CONall.pval,'BHFDR',true);
    CONall = CONall(1:100,:);
    GOtableCON = table();
    GOtableCON.GOcategory = CONall.GOID;
    GOtableCON.Description = CONall.GOName;
    GOtableCON.NumGenes = CONall.numGenes;
    num_dig = 15;
    GOtableCON.Pval = round(CONall.pval*(10^num_dig))/(10^num_dig);
    GOtableCON.Pval_corr = round(CONall.corr_pval*(10^num_dig))/(10^num_dig);
    
    % save file as .csv
    fileOUTname = sprintf('enrichment/output/filter/ermineJresults_%s.csv', type);
    writetable(GOtableCON,fileOUTname)
end

end




