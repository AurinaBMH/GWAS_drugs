function LD_SNPs = SQL_SNP_LD_SNP(mySNP,threshold)
% Retrieves a set of SNPs that are LD to a given input SNP
%-------------------------------------------------------------------------------

if nargin < 2
    threshold = 0;
end

% Connect:
dbc = SQL_ConnectToLD();

if threshold==0
    % Put the query:
    selectText = sprintf(['SELECT SNP2 FROM LD_rel WHERE SNP1=''%s'''],mySNP);
    [LD_SNPs,b,~,~] = mysql_dbquery(dbc,selectText);
else
    % Put the query:
    selectText = sprintf(['SELECT SNP2,r2 FROM LD_rel WHERE SNP1=''%s'''],mySNP);
    [SNPr2,b,~,~] = mysql_dbquery(dbc,selectText);
    r2 = [SNPr2{:,2}];
    strongEnough = (r2 > threshold);
    % SNP_list = SNPr2(:,1);
    LD_SNPs = SNPr2(strongEnough,1);
end

% Check if empty (no matches):
if isempty(LD_SNPs)
    warning('No matching SNPs found for %u',mySNP)
end

SQL_closedatabase(dbc);

end
