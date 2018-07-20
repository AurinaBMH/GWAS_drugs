function geneList = SQL_genesForSNPs(SNPlist)
% Retrieves a list of genes matching a given set of input SNPs
% Taken from the SNPtogene table in the geneLD mySQL database...
%-------------------------------------------------------------------------------

% Connect:
dbc = SQL_ConnectToLD();

% Put the query:
if iscell(SNPlist)
    % Write comma-delimited set of SNPids
    SNP_list_string = BF_cat(SNPlist,',','''');
    selectText = sprintf(['SELECT geneName FROM SNPtogene WHERE SNPname IN (%s)'],SNP_list_string);
else
    selectText = sprintf(['SELECT geneName FROM SNPtogene WHERE SNPname = ''%s'''],SNPlist);
end

[geneList,~,~,emsg] = mysql_dbquery(dbc,selectText);

if ~isempty(emsg)
    keyboard
end

% Check if empty (no matches):
if isempty(geneList)
    geneList = {};
    warning('No genes found for any of the %u input SNPs',length(SNPlist))
else
    geneList = unique(geneList);
end

if length(geneList)==1
    geneList = geneList{1};
end

% Close the database connection:
SQL_closedatabase(dbc);

end
