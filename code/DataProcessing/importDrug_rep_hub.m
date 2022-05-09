function data = importDrug_rep_hub(filename, dataLines)
%IMPORTFILE Import data from a text file
%  DRUGREPURPOSINGHUBDATABASE20200730 = IMPORTFILE(FILENAME) reads data
%  from text file FILENAME for the default selection.  Returns the data
%  as a table.
%
%  DRUGREPURPOSINGHUBDATABASE20200730 = IMPORTFILE(FILE, DATALINES)
%  reads data for the specified row interval(s) of text file FILENAME.
%  Specify DATALINES as a positive scalar integer or a N-by-2 array of
%  positive scalar integers for dis-contiguous row intervals.
%
%  Example:
%  Drugrepurposinghubdatabase20200730 = importfile("/Users/aurinaa/Google_drive/PostDoc/projects/GWASdrugs/data/TREATMENTlists/Drug_repurposing_hub_database/Drug_repurposing_hub_database_20200730.txt", [11, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 04-Aug-2020 11:37:22

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [11, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["pert_iname", "clinical_phase", "moa", "target", "disease_area", "indication"];
opts.VariableTypes = ["string", "categorical", "string", "string", "categorical", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["pert_iname", "moa", "target", "indication"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["pert_iname", "clinical_phase", "moa", "target", "disease_area", "indication"], "EmptyFieldRule", "auto");

% Import the data
data = readtable(filename, opts);

end