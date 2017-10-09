function [AdjPPI,geneNames,PPIN] = PPINImport(evidenceThreshold)
% Import PPIN data from STRING
%-------------------------------------------------------------------------------

if nargin < 1
    evidenceThreshold = 0;
end

% Set where to save to:
fileNameSave = sprintf('PPIN_processed_th%u.mat',evidenceThreshold);
fileNameSave = fullfile('Data',fileNameSave);

% Read in data:
fileName = '6_PPIN_STRINGv10.5.csv';
fprintf(1,'Reading in data from %s...',fileName);
fid = fopen(fileName,'r');
C = textscan(fid,'%s%s%u','Delimiter',',','HeaderLines',1);
fclose(fid);
gene1 = C{1};
gene2 = C{2};
evidenceScore = C{3};
fprintf(1,' Done for %u interactions :-O',length(gene1));

isGood = cellfun(@(x)~isempty(x),gene1) & cellfun(@(x)~isempty(x),gene2);
highEvidence = (evidenceScore>evidenceThreshold);
keepEdge = (isGood & highEvidence);
numInteractions = sum(keepEdge);
fprintf(1,'Filtered on evidence -> %u edges in the PPIN\n',numInteractions);

PPIN = [gene1(keepEdge),gene2(keepEdge)];
clear('evidenceScore','gene1','gene2');

save(fileNameSave,'PPIN');

%-------------------------------------------------------------------------------
% All genes in the PPIN:
%-------------------------------------------------------------------------------
fprintf(1,'Determining unique genes...\n');
geneNames = unique(vertcat(PPIN(:,1),PPIN(:,2)));
numGenes = length(geneNames);
fprintf(1,'%u unique genes in the PPIN\n',numGenes);

% Save:
save(fileNameSave,'geneNames','-append');
fprintf(1,'(saved to %s)\n',fileNameSave);

% Convert gene names to indices:
fprintf(1,['Constructing sparse symmetric adjacency matrix (%ux%u);',...
            ' %u edges for comparison...\n'],numGenes,numGenes,numInteractions);
ii = zeros(numInteractions,1);
jj = zeros(numInteractions,1);
for k = 1:numInteractions
    ii(k) = strcmp(geneNames,PPIN{k,1});
    jj(k) = strcmp(geneNames,PPIN{k,2});
    % AdjPPI(ii,jj) = 1;
end
fprintf(1,'Generating a sparse matrix from the gene-matched indices:\n');
AdjPPI = sparse(ii,jj,true,numGenes,numGenes);
fprintf(1,'Symmetrizing the sparse matrix...');
AdjPPI = (AdjPPI | AdjPPI');
fprintf(1,' Done.\n');

%-------------------------------------------------------------------------------
% Save to .mat file:
save(fileNameSave,'AdjPPI','-append');

end
