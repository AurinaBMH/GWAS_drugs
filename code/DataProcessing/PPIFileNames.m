function fileNames = PPIFileNames(doWeighted,evidenceThreshold,whatInput)
% Turns a set of processing settings into a set of relevant filenames
%-------------------------------------------------------------------------------

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
switch whatInput
case 'HGNCmatch'
    preText = 'PPI_HGNC';
case 'proteins'
    preText = 'PPI_protein';
case 'original'
    preText = 'PPI_gene';
otherwise
    error('Unknown PPI data source: ''%s''',whatInput);
end

if doWeighted
    extraText = '_w';
else
    extraText = sprintf('_th%u',evidenceThreshold);
end

fileNames = cell(4,1);
fileNames{1} = sprintf('%s_processed%s.mat',preText,extraText);
fileNames{2} = sprintf('%s_geneLabels%s.mat',preText,extraText);
fileNames{3} = sprintf('%s_Adj%s.mat',preText,extraText);
fileNames{4} = sprintf('%s_Dist%s.mat',preText,extraText);

% Put them all in the DataOutput directory:
for i = 1:4
    fileNames{i} = fullfile('DataOutput_2024',fileNames{i});
end

end
