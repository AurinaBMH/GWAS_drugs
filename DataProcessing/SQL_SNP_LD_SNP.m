function LD_SNPs = SQL_SNP_LD_SNP(mySNP,threshold,lookTwice)
% Retrieves a set of SNPs that are LD to a given input SNP
%-------------------------------------------------------------------------------

if nargin < 2
    threshold = 0;
end
if nargin < 3
    lookTwice = false;
end

% Connect:
dbc = SQL_ConnectToLD();

if threshold==0
    % Put the query:
    selectText = sprintf('SELECT SNP2 FROM LD_rel WHERE SNP1=''%s''',mySNP);
    [LD_SNPs,~,~,emsg] = mysql_dbquery(dbc,selectText);
    % The other way around:
    if lookTwice
        selectText = sprintf('SELECT SNP1 FROM LD_rel WHERE SNP2=''%s''',mySNP);
        [LD_SNPs_2,~,~,emsg] = mysql_dbquery(dbc,selectText);
        LD_SNPs = union(LD_SNPs,LD_SNPs_2);
    end
else
    % Put the query:
    selectText = sprintf('SELECT SNP2,r2 FROM LD_rel WHERE SNP1=''%s''',mySNP);
    [SNPr2,~,~,emsg] = mysql_dbquery(dbc,selectText);

    if lookTwice
        selectText = sprintf('SELECT SNP1,r2 FROM LD_rel WHERE SNP2=''%s''',mySNP);
        [SNPr2_2,~,~,emsg] = mysql_dbquery(dbc,selectText);
        SNPr2 = vertcat(SNPr2,SNPr2_2);
    end

    if isempty(SNPr2)
        LD_SNPs = {};
    else
        % Consider unique LD SNPs:
        [~,ia] = unique(SNPr2(:,1));
        SNPr2 = SNPr2(ia,:);
        % Take r2 values and threshold:
        r2 = [SNPr2{:,2}];
        strongEnough = (r2 > threshold);
        % SNP_list = SNPr2(:,1);
        LD_SNPs = SNPr2(strongEnough,1);
    end
end

if ~isempty(emsg)
    keyboard
end

% Check if empty (no matches):
if isempty(LD_SNPs)
    LD_SNPs = {};
    % warning('No other SNPs are within-threshold LD with %s',mySNP)
end

SQL_closedatabase(dbc);

end
