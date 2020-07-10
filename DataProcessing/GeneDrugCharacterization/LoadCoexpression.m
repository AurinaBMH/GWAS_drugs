function [geneCoexp,probeInformation] = LoadCoexpression(zscoreOrSigmoid)
% Loads in cortical coexpression data for genes
% Data kindly processed and provided by Aurina Arnatkeviciute
%-------------------------------------------------------------------------------

if nargin < 1
    zscoreOrSigmoid = 'sigmoid';
    fprintf(1,'Sigmoid normalization by default\n');
end
%-------------------------------------------------------------------------------

switch zscoreOrSigmoid
case 'zscore'

case 'sigmoid'
    % this file is the renamed version of
    % 100DS360scaledRobustSigmoidNSGDSQC1Lcortex_ROI_NOdistCorrSurface.mat
    fileName = '360geneExpression.mat';
end
fprintf(1,'Loading data from %s\n',fileName);
load(fileName,'parcelExpression','probeInformation');
fprintf(1,'Loaded; now calculating gene x gene coexpression\n ...');
geneCoexp = corr(parcelExpression(:,2:end), 'type', 'Spearman', 'rows','complete');
fprintf(1,' Done.\n');

end
