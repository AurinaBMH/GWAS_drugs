function params = SetDefaultParams()

%-------------------------------------------------------------------------------
% PPI network:
%-------------------------------------------------------------------------------
% Whether to use evidence scores as a proxy weight for network:
params.doWeighted = true;

params.whatNorm = 2; 

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
params.geneScore = 'P';

%-------------------------------------------------------------------------------
% Scoring similarity between two sets of SNPs/genes
%-------------------------------------------------------------------------------
params.whatScore = 'weightedSum'; %'Kendall', 'weightedSum'

%-------------------------------------------------------------------------------
% How to select genes for PPI mapping
%-------------------------------------------------------------------------------
params.whatThreshold = 'BF';
params.whatDrugTargets = '2020';

params.whatTargets = 'active';

params.whatGWAS = {'ADHD', 'MDD2', 'SCZ', 'BIP2', 'DIABETES', 'HF', 'IBD', 'RA'};
params.whatDiseases_Treatment = {'ADHD', 'BIP', 'SCZ', 'MDD', 'DIABETES', 'IBD', 'HF', 'RA', 'gastro', 'pulmonary'};

params.whatANNOT = {'MAGMAdefault', 'Adult_brain', 'Fetal_brain', 'Neuro', 'Astro', ...
    'eQTLbrain', 'eQTLWhole_Blood', 'eQTLLiver', 'eQTLHeart_Left_Ventricle', 'eQTLPancreas', ...
    'eQTLSmall_Intestine_Terminal_Ileum', 'eQTLColon_Transverse' 'eQTLColon_Sigmoid', 'eQTLAdipose_Subcutaneous' 'eQTLAdipose_Visceral_Omentum'};
end
