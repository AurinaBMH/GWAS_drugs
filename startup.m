% Make folders that will be used to save data
cd ..

mkdir DataOutput_2022
mkdir enrichment_2021
mkdir enrichment_2022
mkdir figures_2021
mkdir figures_2022

% Add paths required for the project (ignoring hidden, including version control)
files = dir;
directories = files([files.isdir]);
directories(strmatch('.',{files([files.isdir]).name})) = []; % remove hidden
paths = arrayfun(@(x)fullfile(directories(x).folder,directories(x).name),1:length(directories),'UniformOutput',false);
for j = 1:length(paths)
    addpath(genpath(paths{j}))
end

%-------------------------------------------------------------------------------
% load('processedData.mat');
%-------------------------------------------------------------------------------
% fprintf(1,'Adding path to MatlabmySQL toolbox\n');
% addpath('/Users/benfulcher/DropboxSydneyUni/CodeToolboxes/MatlabmySQL/')
