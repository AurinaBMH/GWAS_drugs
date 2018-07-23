function [AdjPPI,geneNames,PPIN] = PPINImport(evidenceThreshold,doWeighted)
% Import PPIN data from STRING
%-------------------------------------------------------------------------------

if nargin < 1
    evidenceThreshold = 0.0;
end
if nargin < 2
    doWeighted = true;
end

%-------------------------------------------------------------------------------
% Set where to save to:
% Save outputs across 3 different files:
if doWeighted
    extraText = '_w';
else
    extraText = sprintf('_th0%u',evidenceThreshold*10);
end
fileNameSave1 = sprintf('PPIN_processed%s.mat',extraText);
fileNameSave1 = fullfile('DataOutput',fileNameSave1);
fileNameSave2 = sprintf('PPI_geneLabels%s.mat',extraText);
fileNameSave2 = fullfile('DataOutput',fileNameSave2);
fileNameSave3 = sprintf('PPI_Adj%s.mat',extraText);
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

if doWeighted
    keepEdge = isGood;
    PPIN = [gene1(keepEdge),gene2(keepEdge)];
    evidenceScore = evidenceScore(keepEdge);
else
    highEvidence = (evidenceScore > evidenceThreshold);
    keepEdge = (isGood & highEvidence);
    PPIN = [gene1(keepEdge),gene2(keepEdge)];
end
numInteractions = sum(keepEdge);
fprintf(1,'Filtered on evidence -> %u edges in the PPIN\n',numInteractions);

clear('gene1','gene2','keepEdge','isGood','highEvidence');

%-------------------------------------------------------------------------------
% All genes in the PPIN:
%-------------------------------------------------------------------------------
fprintf(1,'Determining unique genes...\n');
geneNames = unique(vertcat(PPIN(:,1),PPIN(:,2)));
numGenes = length(geneNames);
fprintf(1,'%u unique genes in the PPIN\n',numGenes);

% Save:
save(fileNameSave2,'geneNames');
fprintf(1,'(saved to %s)\n',fileNameSave2);

% Convert gene names to indices:
fprintf(1,['Constructing sparse list of edges (%ux%u);',...
            ' as indices across %u edges...\n'],numGenes,numGenes,numInteractions);
indx1 = zeros(numInteractions,1,'uint8'); % save memory/time by preallocating as 8-bit positive integers
indx2 = zeros(numInteractions,1,'uint8'); % save memory/time by preallocating as 8-bit positive integers
parfor k = 1:numInteractions
    % there can only be one match (since geneNames are unique)
    indx1(k) = find(strcmp(geneNames,PPIN{k,1}),1);
    indx2(k) = find(strcmp(geneNames,PPIN{k,2}),1);
end
save(fileNameSave3,'indx1','indx2');
fprintf(1,'Indices saved to %s\n',fileNameSave3);
clear('PPIN');
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Create the sparse matrix of the interactions:
%-------------------------------------------------------------------------------
fprintf(1,'Generating a sparse matrix from the gene-matched indices:\n');
if doWeighted
    AdjPPI = sparse(double(indx1),double(indx2),double(evidenceScore),numGenes,numGenes);
else
    AdjPPI = sparse(double(indx1),double(indx2),1,numGenes,numGenes);
    fprintf(1,'Symmetrizing the sparse matrix...');
    AdjPPI = (AdjPPI | AdjPPI');
    fprintf(1,' Done.\n');
end
save(fileNameSave3,'AdjPPI','-append');
fprintf(1,'Saved (sparse) Adj to %s\n',fileNameSave3);

end
