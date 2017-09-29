function [AdjPPI,geneNames,PPIN] = PPINImport(evidenceThreshold)
% Import PPIN data from STRING
%-------------------------------------------------------------------------------

if nargin < 1
    evidenceThreshold = 0;
end

fid = fopen('6_PPIN_STRINGv10.5.csv','r');
C = textscan(fid,'%s%s%u','Delimiter',',','HeaderLines',1);
fclose(fid);
gene1 = C{1};
gene2 = C{2};
evidenceScore = C{3};

isGood = cellfun(@(x)~isempty(x),gene1) & cellfun(@(x)~isempty(x),gene2);
highEvidence = (evidenceScore>evidenceThreshold);
keepEdge = (isGood & highEvidence);
numInteractions = sum(keepEdge);
fprintf(1,'%u PPIN edges\n',numInteractions);

PPIN = [gene1(keepEdge),gene2(keepEdge)];

% All genes in the PPIN:
geneNames = unique([unique(PPIN(:,1));unique(PPIN(:,2))]);
numGenes = length(geneNames);

AdjPPI = sparse(numGenes,numGenes);
for k = 1:numInteractions
    ii = strcmp(geneNames,PPIN{k,1});
    jj = strcmp(geneNames,PPIN{k,2});
    AdjPPI(ii,jj) = 1;
end
AdjPPI = (AdjPPI | AdjPPI');

%-------------------------------------------------------------------------------
% Save to .mat file:
fileName = sprintf('PPIN_processed_th%u.mat',evidenceThreshold);
save(fullfile('Data',fileName),'AdjPPI','geneNames','PPIN');

end
