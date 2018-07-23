function distMatrix = ComputePPIDist(evidenceThreshold,doWeighted)
% Compute pairwise distances on the PPI interaction network (from STRING)
%-------------------------------------------------------------------------------
if nargin < 1
    evidenceThreshold = 0;
end
if nargin < 2
    doWeighted = true;
end

%-------------------------------------------------------------------------------
% Load in PPI data:
if doWeighted
    extraText = '_w';
else
    extraText = sprintf('_th0%u',evidenceThreshold*10);
end
fileNameLoad = sprintf('PPI_Adj%s.mat',extraText);
fileNameSave = sprintf('PPI_Dist%s.mat',extraText);
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
save(fullfile('Data',fileNameSave),'distMatrix');
fprintf(1,'Saved pairwise network distance data to %s\n',fileNameSave);

end
