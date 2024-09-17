function treatmentlistBIPlithium = import_drug_class_file(filename, dataLines)
%IMPORTFILE Import data from a text file
%  TREATMENTLISTBIPLITHIUM = IMPORTFILE(FILENAME) reads data from text
%  file FILENAME for the default selection.  Returns the data as a cell
%  array.
%
%  TREATMENTLISTBIPLITHIUM = IMPORTFILE(FILE, DATALINES) reads data for
%  the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  treatmentlistBIPlithium = importfile("/Users/aurinaa/Documents/PostDoc/projects/GWASdrugs/data/TREATMENTlists/Drug_list_disorders/treatment_list_BIP_lithium.txt", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 17-Sep-2024 11:26:46

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 1);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "";

% Specify column names and types
opts.VariableNames = "lithiumCarbonate";
opts.VariableTypes = "char";

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";

% Specify variable properties
opts = setvaropts(opts, "lithiumCarbonate", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "lithiumCarbonate", "EmptyFieldRule", "auto");

% Import the data
treatmentlistBIPlithium = readtable(filename, opts);

%% Convert to output type
treatmentlistBIPlithium = table2cell(treatmentlistBIPlithium);
numIdx = cellfun(@(x) ~isnan(str2double(x)), treatmentlistBIPlithium);
treatmentlistBIPlithium(numIdx) = cellfun(@(x) {str2double(x)}, treatmentlistBIPlithium(numIdx));
end