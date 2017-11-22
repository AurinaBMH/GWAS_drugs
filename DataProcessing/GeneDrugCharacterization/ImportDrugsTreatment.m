function drugTable = ImportDrugsTreatment()
% Janette has prepared a data file containing drugs currently used to treat
% a range of disorders
%-------------------------------------------------------------------------------

%% Import the data
[~, ~, raw] = xlsread('Drug_identifier.xlsx','Drug_identifier');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
stringVectors = string(raw(:,[1,2,3,4,5,6]));
stringVectors(ismissing(stringVectors)) = '';
raw = raw(:,[7,8,9,10,11,12,13,14,15,16,17]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Create table
drugTable = table;

%% Allocate imported array to column variable names
drugTable.Protein_DrugIDs = stringVectors(:,1);
drugTable.DrugIDs_all = stringVectors(:,2);
drugTable.Indication_MediHPS = stringVectors(:,3);
drugTable.DrugIDs_unnest = stringVectors(:,4);
drugTable.Drug_CommonName_DB = stringVectors(:,5);
drugTable.Action = categorical(stringVectors(:,6));
drugTable.N01_Anesthetics = logical(data(:,1));
drugTable.N02_Analgesics = logical(data(:,2));
drugTable.N03_antispileptics = logical(data(:,3));
drugTable.N04_antiparkinson = logical(data(:,4));
drugTable.N05A_Antipsychotics = logical(data(:,5));
drugTable.N05B_Antianxiety = logical(data(:,6));
drugTable.N05C_Hypotonicssedative = logical(data(:,7));
drugTable.N06A_Antidepressants = logical(data(:,8));
drugTable.N06B_C_Psychostimulants = logical(data(:,9));
drugTable.N06D_antidementia = logical(data(:,10));
drugTable.N07_otherCNSdrugs = logical(data(:,11));

%-------------------------------------------------------------------------------
% Add more:
%---schiz
isSchiz = strcmp(drugTable.Indication_MediHPS,'Schizophrenic disorders');
drugTable.isSchiz = isSchiz;
%---ADHD:
isADHD = startsWith(drugTable.Indication_MediHPS,'Attention');
drugTable.isADHD = isADHD;

end
