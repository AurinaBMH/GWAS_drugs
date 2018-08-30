function [AdjPPI,geneNames,PPIN] = PPINImport(doWeighted,evidenceThreshold,whatInput)
% Import PPIN data from STRING
%-------------------------------------------------------------------------------

% Check inputs:
if nargin < 1
    doWeighted = false;
end
if nargin < 2
    evidenceThreshold = 400;
end
if nargin < 3
    whatInput = 'HGNCmatch';
end

%-------------------------------------------------------------------------------
% Filenames to save processed data to:
fileNameSave = PPIFileNames(doWeighted,evidenceThreshold,whatInput);

%-------------------------------------------------------------------------------
% Read in PPI data:
switch whatInput
case 'HGNCmatch'
    % Latest matching to HGNC genes by Janette (July 2018):
    fileName = 'PPI_conversionToGenes.csv';
case 'proteins'
    % Mapping to proteins:
    fileName = '9606.protein.links.v10.5.txt';
case 'original'
    % Original STRING nomenclature:
    fileName = '6_PPIN_STRINGv10.5.csv';
end
fprintf(1,'Constructing protein-protein interactions using %s input: ''%s''\n',...
                whatInput,fileName);
fid = fopen(fileName,'r');
C = textscan(fid,'%s%s%u','Delimiter',',','HeaderLines',1);
fclose(fid);
gene1 = C{1};
gene2 = C{2};
evidenceScore = C{3};
clear('C');
fprintf(1,'Read in %u interactions :-O\n',length(gene1));

% Strip whitespace:
gene1 = strtrim(gene1);
gene2 = strtrim(gene2);

%-------------------------------------------------------------------------------
% Construct a PPI network, PPIN:
isGood = cellfun(@(x)~isempty(x),gene1) & cellfun(@(x)~isempty(x),gene2);
if strcmp(whatInput,'HGNCmatch')
    % Remove 'uncharacterised' entries:
    isUncharacterised = strcmp(gene1,'uncharacterised') | strcmp(gene2,'uncharacterised');
    isGood(isUncharacterised) = false;
    fprintf(1,'%u ''uncharacterised'' entries filtered out\n',sum(isUncharacterised));
end
if doWeighted
    keepEdge = isGood;
else
    highEvidence = (evidenceScore > evidenceThreshold);
    keepEdge = (isGood & highEvidence);
end
% Filter edges by keepEdge
PPIN = [gene1(keepEdge),gene2(keepEdge)];
evidenceScore = evidenceScore(keepEdge);
numInteractions = sum(keepEdge);
if doWeighted
    fprintf(1,'Weighted analysis: keeping all %u edges in the PPIN\n',numInteractions);
else
    fprintf(1,'Filtered on evidence threshold %g -> %u edges in the PPIN\n',...
                        evidenceThreshold,numInteractions);
end

clear('gene1','gene2','keepEdge','isGood','highEvidence');

%-------------------------------------------------------------------------------
% Determine the list of genes contained in the PPIN:
fprintf(1,'Determining unique genes...\n');
geneNames = unique(vertcat(PPIN(:,1),PPIN(:,2)));
numGenes = length(geneNames);
fprintf(1,'%u unique genes in the PPIN\n',numGenes);

% Save:
save(fileNameSave{2},'geneNames');
fprintf(1,'(saved to %s)\n',fileNameSave{2});

% Convert gene names to indices:
fprintf(1,['Constructing sparse list of edges (%ux%u);',...
            ' as indices across %u edges...\n'],numGenes,numGenes,numInteractions);
indx1 = zeros(numInteractions,1,'uint16');
indx2 = zeros(numInteractions,1,'uint16');
if length(geneNames) > 2^16
    error('Problem storing indices as unsigned int16 data');
end
parfor k = 1:numInteractions
    % There can only be one match (since geneNames contains unique entries)
    indx1(k) = find(strcmp(PPIN{k,1},geneNames),1);
    indx2(k) = find(strcmp(PPIN{k,2},geneNames),1);
end
save(fileNameSave{3},'indx1','indx2');
fprintf(1,'Indices saved to %s\n',fileNameSave{3});
clear('PPIN');

%-------------------------------------------------------------------------------
% Create a sparse matrix containing all interactions:
fprintf(1,'Generating a sparse matrix from the gene-matched indices:\n');
if doWeighted
    AdjPPI = sparse(double(indx1),double(indx2),double(evidenceScore),numGenes,numGenes);
else
    AdjPPI = sparse(double(indx1),double(indx2),1,numGenes,numGenes);
    fprintf(1,'Symmetrizing the sparse matrix...');
    AdjPPI = (AdjPPI | AdjPPI');
    fprintf(1,' Done.\n');
end
save(fileNameSave{3},'AdjPPI','-append');
fprintf(1,'Saved (sparse) Adj to %s\n',fileNameSave{3});

end
