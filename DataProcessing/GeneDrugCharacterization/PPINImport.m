function PPIN = PPINImport(threshold)
% Import PPIN data from STRING
%-------------------------------------------------------------------------------

if nargin < 1
    threshold = 0;
end

fid = fopen('6_PPIN_STRINGv10.5.csv','r');
C = textscan(fid,'%s%s%u','Delimiter',',','HeaderLines',1);
fclose(fid);
gene1 = C{1};
gene2 = C{2};
evidenceScore = C{3};

isGood = cellfun(@(x)~isempty(x),gene1) & cellfun(@(x)~isempty(x),gene2);
highEvidence = (evidenceScore>threshold);
keepEdge = (isGood & highEvidence);
fprintf(1,'%u PPIN edges\n',sum(keepEdge));

PPIN = [gene1(keepEdge),gene2(keepEdge)];

end
