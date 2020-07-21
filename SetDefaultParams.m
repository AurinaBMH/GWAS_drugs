function params = SetDefaultParams()

%-------------------------------------------------------------------------------
% PPI network:
%-------------------------------------------------------------------------------
% Whether to use evidence scores as a proxy weight for network:
params.doWeighted = true;

% Evidence threshold for including PPI interactions in binary analyses:
params.PPINevidenceThreshold = 0;

%-------------------------------------------------------------------------------
% SNP-gene mapping
%-------------------------------------------------------------------------------
% Threshold for linking SNPs through LD relationships:
% params.LDthreshold = 0.5;

%-------------------------------------------------------------------------------
% Score for genes - could be ZSTAT, could be P
%-------------------------------------------------------------------------------
params.geneScore = 'ZSTAT'; 

%-------------------------------------------------------------------------------
% Scoring similarity between two sets of SNPs/genes
%-------------------------------------------------------------------------------
params.whatScore = 'weightedSum'; %'Kendall', 'weightedSum'

%-------------------------------------------------------------------------------
% How to select genes for PPI mapping
%-------------------------------------------------------------------------------
params.whatThreshold = 'BF'; 
end
