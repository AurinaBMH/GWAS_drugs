function [AdjPPI,geneNames,PPIN] = PPINImport(evidenceThreshold)
% Import PPIN data from STRING
%-------------------------------------------------------------------------------

if nargin < 1
    evidenceThreshold = 0.0;
end

%-------------------------------------------------------------------------------
% Set where to save to:
% Save outputs across 3 different files:
fileNameSave1 = sprintf('PPIN_processed_th0%u.mat',evidenceThreshold*10);
fileNameSave1 = fullfile('DataOutput',fileNameSave1);
fileNameSave2 = sprintf('PPI_geneLabels_th0%u.mat',evidenceThreshold*10);
fileNameSave2 = fullfile('DataOutput',fileNameSave2);
fileNameSave3 = sprintf('PPI_Adj_th0%u.mat',evidenceThreshold*10);
fileNameSave3 = fullfile('DataOutput',fileNameSave3);

%-------------------------------------------------------------------------------
% Read in data:
fileName = '6_PPIN_STRINGv10.5.csv';
fprintf(1,'Reading in data from %s...',fileName);
fid = fopen(fileName,'r');
C = textscan(fid,'%s%s%u','Delimiter',',','HeaderLines',1);
fclose(fid);
gene1 = C{1};
gene2 = C{2};
evidenceScore = C{3};
clear('C');
fprintf(1,' Done for %u interactions :-O\n',length(gene1));

isGood = cellfun(@(x)~isempty(x),gene1) & cellfun(@(x)~isempty(x),gene2);
highEvidence = (evidenceScore>evidenceThreshold);
keepEdge = (isGood & highEvidence);
numInteractions = sum(keepEdge);
fprintf(1,'Filtered on evidence -> %u edges in the PPIN\n',numInteractions);

PPIN = [gene1(keepEdge),gene2(keepEdge)];
clear('evidenceScore','gene1','gene2','keepEdge','isGood','highEvidence');

% This may be too big for a .mat file...?
% fprintf(1,'Saving filtered data to %s...',fileNameSave1);
% save(fileNameSave1,'PPIN','-v7.3');
% fprintf(1,' Saved.\n',fileNameSave1);

%-------------------------------------------------------------------------------
% All genes in the PPIN:
%-------------------------------------------------------------------------------
fprintf(1,'Determining unique genes...\n');
geneNames = unique(vertcat(PPIN(:,1),PPIN(:,2)));
numGenes = length(geneNames);
fprintf(1,'%u unique genes in the PPIN\n',numGenes);

% Save:
save(fileNameSave2,'geneNames','-v7.3');
fprintf(1,'(saved to %s)\n',fileNameSave2);

% Convert gene names to indices:
fprintf(1,['Constructing sparse symmetric adjacency matrix (%ux%u);',...
            ' %u edges for comparison...\n'],numGenes,numGenes,numInteractions);
indx1 = zeros(numInteractions,1,'uint8'); % save memory/time by preallocating as 8-bit positive integers
indx2 = zeros(numInteractions,1,'uint8'); % save memory/time by preallocating as 8-bit positive integers
for k = 1:numInteractions
    % there can only be one match (since geneNames are unique)
    indx1(k) = find(strcmp(geneNames,PPIN{k,1}),1);
    indx2(k) = find(strcmp(geneNames,PPIN{k,2}),1);
end
save(fileNameSave3,'indx1','indx2','-v7.3');
fprintf(1,'Indices saved to %s\n',fileNameSave3);
clear('PPIN');
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Create the sparse matrix of the interactions:
%-------------------------------------------------------------------------------
fprintf(1,'Generating a sparse matrix from the gene-matched indices:\n');
AdjPPI = sparse(double(indx1),double(indx2),1,numGenes,numGenes);
fprintf(1,'Symmetrizing the sparse matrix...');
AdjPPI = (AdjPPI | AdjPPI');
fprintf(1,' Done.\n');
save(fileNameSave3,'AdjPPI','-append');
fprintf(1,'Saved symmetrized (sparse) Adj to %s\n',fileNameSave3);

end
