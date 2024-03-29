% Make folders that will be used to save data
mkdir data
mkdir DataOutput_2024
mkdir figures_2021
mkdir figures_2024
mkdir enrichment_2024
mkdir enrichment_2024/output
mkdir enrichment_2021
mkdir enrichment_2021/output

% Add paths required for the project (ignoring hidden, including version control)
files = dir;
directories = files([files.isdir]);
directories(strmatch('.',{files([files.isdir]).name})) = []; % remove hidden
paths = arrayfun(@(x)fullfile(directories(x).folder,directories(x).name),1:length(directories),'UniformOutput',false);
for j = 1:length(paths)
    addpath(genpath(paths{j}))
end

