function distMatrix = ComputePPIDist(evidenceThreshold)
% Import PPIN data from STRING
%-------------------------------------------------------------------------------
if nargin < 1
    evidenceThreshold = 0;
end
fileNameLoad = sprintf('PPI_Adj_th%u.mat',evidenceThreshold);
fileNameSave = sprintf('PPI_Dist_th%u.mat',evidenceThreshold);
load(fileNameLoad,'AdjPPI')
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Compute pairwise distances:
fprintf(1,'Computing pairwise distances between %u genes...\n',size(AdjPPI,1));
distMatrix = distance_bin(AdjPPI);

%-------------------------------------------------------------------------------
% Save out:
save(fullfile('Data',fileNameSave),'distMatrix');
fprintf(1,'Saved pairwise network distance data to %s\n',fileNameSave);

end
