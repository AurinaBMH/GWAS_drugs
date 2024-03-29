function BioMartRESULTS = importBIOMARTfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  BIOMARTRESULTS = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  BIOMARTRESULTS = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  BioMartRESULTS = importfile("/Users/aurinaa/Google_drive/PostDoc/projects/GWASdrugs/rawData/BioMart_RESULTS.txt", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 25-Jun-2020 11:58:20

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 9);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["GenestableID", "GenestableIDversion", "TranscriptstableID", "TranscriptstableIDversion", "ProteinstableID", "ProteinstableIDversion", "Genename", "HGNCsymbol", "HGNCID"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["GenestableID", "GenestableIDversion", "TranscriptstableID", "TranscriptstableIDversion", "ProteinstableID", "ProteinstableIDversion", "Genename", "HGNCsymbol", "HGNCID"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["GenestableID", "GenestableIDversion", "TranscriptstableID", "TranscriptstableIDversion", "ProteinstableID", "ProteinstableIDversion", "Genename", "HGNCsymbol", "HGNCID"], "EmptyFieldRule", "auto");

% Import the data
BioMartRESULTS = readtable(filename, opts);

end