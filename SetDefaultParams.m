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
% AHBA measure
%-------------------------------------------------------------------------------
params.AHBAmeasure = 'z'; % could be -z, abs(z)

%-------------------------------------------------------------------------------
% Score for genes - could be ZSTAT, could be P
%-------------------------------------------------------------------------------
params.geneScore = 'P';
params.DOthreshold = false; 

%-------------------------------------------------------------------------------
% Scoring similarity between two sets of SNPs/genes
%-------------------------------------------------------------------------------
params.whatScore = 'weightedSum'; %'Kendall', 'weightedSum'

%-------------------------------------------------------------------------------
% How to select genes for PPI mapping
%-------------------------------------------------------------------------------
params.whatThreshold = 'BF';
params.whatDrugTargets = '2020';
params.whatTargets = 'all';
params.numNull = 5000; 

% keep only psychiatric drugs+ diabetes
params.whatGWAS = {'ADHD', 'MDD2', 'SCZ', 'BIP2', 'DIABETES', 'HF', 'IBD', 'RA'};

params.whatDiseases_Treatment = {'ADHD', 'BIP', 'SCZ', 'MDD', 'DIABETES'};
params.whatDiseases_Treatment_label = {'ADHD', 'Bipolar disorder', 'Schizophrenia', 'Major depression', 'Diabetes'};

params.whatDiseases_Treatment_body = {'DIABETES', 'IBD', 'HF', 'RA', 'gastro', 'pulmonary'};
params.whatDiseases_Treatment_label_body = {'Diabetes', 'IBD', 'Heart failure', 'Rheumatoid arthritis', 'Gastroentherology', 'Pulmonology'};

params.whatDiseases_Treatment_ALL = unique([params.whatDiseases_Treatment, params.whatDiseases_Treatment_body], 'stable');
params.whatDiseases_Treatment_label_ALL = unique([params.whatDiseases_Treatment_label, params.whatDiseases_Treatment_label_body], 'stable');

params.whatANNOT_psych = {'MAGMAdefault', 'Adult_brain', 'Fetal_brain', 'Neuro', 'Astro', ...
    'eQTLbrain', 'eQTLWhole_Blood', 'eQTLPancreas', 'eQTLLiver', ...
    'AllenMapped', 'AlleneQTLbrain', ...
    'PPI_eQTLbrain_th0', 'PPI_eQTLbrain_th400', 'PPI_eQTLbrain_th600', 'PPI_eQTLbrain_th900', ...
    'PPI_mapped_th0', 'PPI_mapped_th400', 'PPI_mapped_th600', 'PPI_mapped_th900'};

params.whatANNOT_body = {'MAGMAdefault', 'eQTLWhole_Blood', 'eQTLLiver', 'eQTLHeart_Left_Ventricle', 'eQTLPancreas', ...
    'eQTLSmall_Intestine_Terminal_Ileum', 'eQTLColon_Transverse' 'eQTLColon_Sigmoid', 'eQTLAdipose_Subcutaneous' 'eQTLAdipose_Visceral_Omentum', ...
    'AllenMapped', 'AlleneQTLbrain', ...
    'PPI_eQTLbrain_th0', 'PPI_eQTLbrain_th400', 'PPI_eQTLbrain_th600', 'PPI_eQTLbrain_th900', ...
    'PPI_mapped_th0', 'PPI_mapped_th400', 'PPI_mapped_th600', 'PPI_mapped_th900'};

params.whatANNOT_all = unique([params.whatANNOT_psych, params.whatANNOT_body]); 

params.whatANNOT_reduced = {'MAGMAdefault', 'PPI_mapped_th600', 'eQTLbrain', 'AlleneQTLbrain'}; 

end
