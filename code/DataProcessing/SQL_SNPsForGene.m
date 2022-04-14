function SNP_list = SQL_SNPsForGene(geneName)
% Retrieves a list of SNPs for a given input gene
% Taken from the SNPtogene table in the geneLD mySQL database...
%-------------------------------------------------------------------------------

% Connect:
dbc = SQL_ConnectToLD();

% Put the query:
selectText = sprintf(['SELECT SNPname FROM SNPtogene WHERE geneName=''%s'''],geneName);
[SNP_list,b,~,~] = mysql_dbquery(dbc,selectText);

% Check if empty (no matches):
if isempty(SNP_list)
    warning('No match found for %u',geneName)
end

% Close the database connection:
SQL_closedatabase(dbc);

end
