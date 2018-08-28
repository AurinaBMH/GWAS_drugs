function distMatrix = ComputePPIDist(evidenceThreshold,doWeighted)
% Compute pairwise distances on the PPI interaction network (from STRING)
%-------------------------------------------------------------------------------
if nargin < 1
    evidenceThreshold = 400;
end
if nargin < 2
    doWeighted = false;
end

%-------------------------------------------------------------------------------
% Load in PPI data:
fileNames = PPIFileNames(doWeighted,evidenceThreshold,'HGNCmatch');
fileNameLoad = fileNames{3};
fileNameSave = fileNames{4};
load(fileNameLoad,'AdjPPI')
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Compute pairwise distances:
numGenes = size(AdjPPI,1);
if doWeighted
    fprintf(1,'Computing weighted pairwise PPI distances between %u genes...\n',numGenes);
    distMatrix = distance_wei(AdjPPI);
else
    fprintf(1,'Computing binary pairwise PPI distances between %u genes...\n',numGenes);
    distMatrix = distance_bin(AdjPPI);
end

%-------------------------------------------------------------------------------
% Save out:
save(fileNameSave,'distMatrix');
fprintf(1,'Saved pairwise network distance data to %s\n',fileNameSave);

end
