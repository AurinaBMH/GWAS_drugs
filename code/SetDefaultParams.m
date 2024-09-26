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
params.whatDrugTargets = '2024';
params.whatTargets = 'all';
params.numNull = 10000; 

% keep only psychiatric drugs+ diabetes
% for 2021 data
params.whatGWAS_2021 = {'ADHD', 'MDD3', 'SCZ', 'BIP2', 'DIABETES', 'HF', 'IBD', 'RA'};
% for 2022 data
%params.whatGWAS_2022 = {'ADHD3','MDD3', 'SCZ3', 'BIP3', 'DIABETES', 'HF', 'IBD', 'RA'};
% for 2024 data
params.whatGWAS_2024 = {'ADHD3','MDD4', 'SCZ3', 'BIP3', 'DIABETES2', 'HF', 'IBD', 'RA'};

params.whatDiseases_Treatment = {'ADHD', 'BIP', 'SCZ', 'MDD', 'DIABETES'};
params.whatDiseases_Treatment_label = {'ADHD', 'Bipolar disorder', 'Schizophrenia', 'Major depression', 'Diabetes'};

params.whatDiseases_Treatment_body = {'DIABETES', 'IBD', 'HF', 'RA', 'gastro', 'pulmonary'};
params.whatDiseases_Treatment_label_body = {'Diabetes', 'IBD', 'Heart failure', 'Rheumatoid arthritis', 'Gastroentherology', 'Pulmonology'};

params.whatDiseases_Treatment_ALL = unique([params.whatDiseases_Treatment, params.whatDiseases_Treatment_body], 'stable');
params.whatDiseases_Treatment_label_ALL = unique([params.whatDiseases_Treatment_label, params.whatDiseases_Treatment_label_body], 'stable');

% list treatment classes BIP
params.whatDiseases_Treatment_classes = {'BIP_lithium', 'BIP_antipsychotic', 'BIP_antidepressant', 'BIP_anticonvulsant'};
params.whatDiseases_Treatment_classes_label = {'BIP lithium', 'BIP antipsychotic', 'BIP antidepressant', 'BIP anticonvulsant'};

% list treatment classes SCZ
params.whatDiseases_Treatment_classes_SCZ = {'SCZ_antidepressant', 'SCZ_atypical_antipsychotic', 'SCZ_typical_antipsychotic'};
params.whatDiseases_Treatment_classes_label_SCZ ={'SCZ antidepressant', 'SCZ atypical antipsychotic', 'SCZ typical antipsychotic'};

% list treatment classes MDD
params.whatDiseases_Treatment_classes_MDD = {'MDD_atypical_antidepressant', 'MDD_atypical_antipsychotic', 'MDD_antipsychotic', 'MDD_MAOI', 'MDD_SNRI', 'MDD_SSRI', 'MDD_TCA'};
params.whatDiseases_Treatment_classes_label_MDD = {'MDD atypical antidepressant', 'MDD atypical antipsychotic', 'MDD antipsychotic', 'MDD MAOI', 'MDD SNRI', 'MDD SSRI', 'MDD TCA'};

% list treatment classes ADHD
params.whatDiseases_Treatment_classes_ADHD = {'ADHD_stimulant', 'ADHD_nonstimulant'};
params.whatDiseases_Treatment_classes_label_ADHD = {'ADHD stimulant', 'ADHD nonstimulant'}; 

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

params.whatANNOT_MAGMA = {'MAGMAdefault', 'Adult_brain', 'Fetal_brain', 'Neuro', 'Astro', ...
    'eQTLbrain', 'eQTLWhole_Blood', 'eQTLLiver', 'eQTLHeart_Left_Ventricle', 'eQTLPancreas', ...
    'eQTLSmall_Intestine_Terminal_Ileum', 'eQTLColon_Transverse' 'eQTLColon_Sigmoid', 'eQTLAdipose_Subcutaneous' 'eQTLAdipose_Visceral_Omentum'};

params.whatANNOT_reduced = {'MAGMAdefault', 'PPI_mapped_th600', 'eQTLbrain', 'AlleneQTLbrain'}; 

params.whatANNOT_reduced5 = {'MAGMAdefault', 'PPI_mapped_th600', 'eQTLbrain', 'AlleneQTLbrain', 'Adult_brain'}; 

end
