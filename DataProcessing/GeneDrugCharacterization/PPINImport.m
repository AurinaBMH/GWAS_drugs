function [AdjPPI,geneNames] = PPINImport(doWeighted,evidenceThreshold,whatInput)
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
        % fileName = 'PPI_conversionToGenes.csv';
        % File generated based on 9606.protein.links.v11.0.txt (June 2020), gene names extracted from biomart using 
        % make_PPI_linkfile.m
        fileName = 'PPIlinks_v11.0.txt'; 
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
 
nPairsR = table2cell(unique(cell2table([gene1,gene2]))); 
fprintf('There are %d entries in the original data, %d pairs are unique\n', length(gene1), length(nPairsR))
clear('nPairsR');

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

[nPairsC,W] = unique(cell2table(PPIN), 'rows','stable'); 
fprintf('There are %d entries in the characterised data, %d pairs are unique\n', size(PPIN,1), size(nPairsC,1))
clear('nPairsC');
% after assigning proteins to genes, some protein pairs (n=36540) become non-unique 
% what to do in that case? 
% take the pair with the highest evidence score; 
DUP = setdiff(1:size(PPIN,1),W); 
isDUP = false(length(evidenceScore),1); 
fprintf('There are %d duplicate entries\n', length(DUP))
fprintf('Select entry with highest evidence score\n')

for d=1:length(DUP)
    V = PPIN(DUP(d),:);
    V1 = find(ismember(PPIN(:,1),V{1})); 
    V2 = find(ismember(PPIN(:,2),V{2})); 
    V12 = intersect(V1,V2); 
    
    % find non-max scores and save their index to exclude in isDUP
    % if all scores are equal, this will select one fo them
    [~, sI] = sort(evidenceScore(V12));
    isDUP(V12(sI(2:end))) = 1;
    fprintf('%d pair processed\n', d); 
end

% remove non-selected duplicates
PPIN(isDUP,:) = []; 
evidenceScore(isDUP) = [];
numInteractions = length(evidenceScore); 

if doWeighted
    fprintf(1,'Weighted analysis: keeping all unique %u edges in the PPIN\n',numInteractions);
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

%AdjPPI = zeros(numGenes, numGenes); 
for k = 1:numInteractions
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
end

    
% for i=1:length(evidenceScore)
%     if doWeighted
% 
%         AdjPPI(double(indx1(i)), double(indx2(i))) = double(evidenceScore(i));
%         
%     else
%         
%         AdjPPI(indx1(i), indx2(i)) = 1;
%         % sparse version generates evidence Scores>1000 for the last gene in the matrix - don't know why;
%         % AdjPPI = sparse(double(indx1),double(indx2),double(evidenceScore),numGenes,numGenes);
%         
%         %AdjPPI = sparse(double(indx1),double(indx2),1,numGenes,numGenes);
%         %fprintf(1,'Symmetrizing the sparse matrix...');
%         %AdjPPI = (AdjPPI | AdjPPI');
%         
%     end
%     
% end
fprintf(1,' Done.\n');

% check if the matrix is symmetric
if issymmetric(AdjPPI)
    fprintf('PPI matrix is symmetric')
else
    fprintf('PPI matrix is not symmetric')
end

save(fileNameSave{3},'AdjPPI','-append');
fprintf(1,'Saved Adj to %s\n',fileNameSave{3});

end
