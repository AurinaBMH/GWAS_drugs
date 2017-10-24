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
    fileName = 'geneXgeneCoexpressionSigmoidSquareform.mat';
end
fprintf(1,'Loading data from %s\n',fileName);
load(fileName,'geneCoexpressionDist','probeInformation');
fprintf(1,'Loaded; now converting to square matrix of Spearman correlation values...');
geneCoexp = 1-squareform(geneCoexpressionDist);
fprintf(1,' Done.\n');

end
