function BIOMARTgeneIDs = importPROTEINfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  BIOMARTGENEIDS = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  BIOMARTGENEIDS = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  BIOMARTgeneIDs = importfile("/Users/aurinaa/Google_drive/PostDoc/projects/GWASdrugs/data/GWASlists/BIOMART_geneIDs.txt", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 08-Jul-2020 13:01:51

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id", "external_gene_name", "hgnc_symbol", "hgnc_id", "entrezgene_id"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id", "external_gene_name", "hgnc_symbol", "hgnc_id"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id", "external_gene_name", "hgnc_symbol", "hgnc_id"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "entrezgene_id", "TrimNonNumeric", true);
opts = setvaropts(opts, "entrezgene_id", "ThousandsSeparator", ",");

% Import the data
BIOMARTgeneIDs = readtable(filename, opts);

end