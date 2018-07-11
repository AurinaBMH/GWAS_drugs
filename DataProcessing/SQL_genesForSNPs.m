function geneList = SQL_genesForSNPs(SNPlist)
% Retrieves a list of genes matching a given set of input SNPs
% Taken from the SNPtogene table in the geneLD mySQL database...
%-------------------------------------------------------------------------------

% Connect:
dbc = SQL_ConnectToLD();

% Put the query:
% Write comma-delimited set of SNPids
SNP_list_string = BF_cat(SNPlist,',','''');
selectText = sprintf(['SELECT geneName FROM SNPtogene WHERE SNPname IN (%s)'],SNP_list_string);
geneList = mysql_dbquery(dbc,selectText);
geneList = unique(geneList);

% Check if empty (no matches):
if isempty(geneList)
    warning('No genes found for any of the %u input SNPs',length(SNPlist))
end

% Close the database connection:
SQL_closedatabase(dbc);

end
