function [indicatorTable,percIndicatorTable, dataTable] = ImportTreatmentLists(normalizeWithinDrugs, whatDrugTargets)
% Import information on gene targets for psychiatric conditions
%-------------------------------------------------------------------------------
% Input parameters:
%-------------------------------------------------------------------------------
if nargin < 1
    normalizeWithinDrugs = true;
end
if nargin < 2
    whatDrugTargets = '2018';
    % 2020 - uses automated AA version
    % 2018 - uses Janett's version from May 2018; 
end

if normalizeWithinDrugs
    fprintf(1,'Weighting genes equally within a drug...\n');
end

%-------------------------------------------------------------------------------
% Whether to treat each drug as equally important, and each gene targeted by a
% drug as equally important for the efficacy of that drug:
whatDiseases = {'ADHD','BIP','SCZ','MDD','pulmonary','cardiology','gastro','diabetes'};
numDiseases = length(whatDiseases);
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%% Initialize variables.
delimiter = ',';
startRow = 1;
endRow = inf;
%% Format for each line of text:
%   column1: text (%q)
%	column2: text (%q)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%[^\n\r]';

dataTable = struct();
switch whatDrugTargets
    case '2018'
        tic
        % use drug targets assigned by Janette in 05/2018
        for k = 1:numDiseases
            whatDisease = whatDiseases{k};
            switch whatDisease
                case 'ADHD'
                    fileName = 'Treatment-list-ADHD-4thMay2018.csv';
                case 'BIP'
                    fileName = 'Treatment-list-BIP-7thMay2018.csv';
                case 'SCZ'
                    fileName = 'Treatment-list-SCZ-4thMay2018.csv';
                case 'MDD'
                    fileName = 'Treatment-list-MDD-7thMay2018.csv';
                case 'pulmonary'
                    fileName = 'Treatment-list-pulmonary-4thMay2018.csv';
                case 'cardiology'
                    fileName = 'Treatment-list-cardiology-4thMay2018.csv';
                case 'gastro'
                    fileName = 'Treatment-list-gastro-4thMay2018.csv';
                case 'diabetes'
                    fileName = 'Treatment-list-diabetes-25thMay2018.csv';
                otherwise
                    error('Unknown disease: ''%s''',whatDisease);
            end
            
            %% Open the text file.
            fileID = fopen(fileName,'r');
            
            %% Read columns of data according to the format.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', 1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            for block=2:length(startRow)
                frewind(fileID);
                dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                for col=1:length(dataArray)
                    dataArray{col} = [dataArray{col};dataArrayBlock{col}];
                end
            end
            fclose(fileID); % Close the text file
            
            % Create structure of tables:
            dataTable.(whatDisease) = table(dataArray{1:end-1},'VariableNames',{'Name','Target'});
            
            numDrugs = length(dataArray{1});
            fprintf(1,'%s has %u drugs\n',whatDisease,numDrugs);
        end
        toc
    case '2020'
        % use drug targets assigned automatically by AA in 08/2020
        % it takes ~10s to run, so load the pre-computed data here
        % dataTable = give_drugTargets();
        load('DataOutput/drugTargets_2020.mat', 'dataTable'); 
end

%-------------------------------------------------------------------------------
% Get a list of all genes mentioned:
geneLists = cell(numDiseases,1);
counts = cell(numDiseases,1);
for k = 1:numDiseases
    whatDisease = whatDiseases{k};
    geneList = dataTable.(whatDisease).Target;
    geneList = cellfun(@(x)strsplit(x,','),geneList,'UniformOutput',false);
    if normalizeWithinDrugs
        numDrugs = length(geneList);
        % Make a weight vector:
        l = cellfun(@length,geneList); % number of genes for each drug
        weightCell = cell(numDrugs,1);
        for j = 1:numDrugs
            weightCell{j} = repmat(1/l(j),1,l(j));
        end
        weights = horzcat(weightCell{:});
    end
    % Genes that appear in drugs for disorder k:
    geneList = horzcat(geneList{:});
    % Remove any whitespace:
    geneList = regexprep(geneList,'\W','');
    
    % Count how many appear:
    frequencyTable = tabulate(geneList);
    geneLists{k} = frequencyTable(:,1);
    numGenesTreat = length(geneLists{k});
    
    % Determine weights on genes:
    if normalizeWithinDrugs
        % The count cell should be weighted by the weights vector:
        counts{k} = zeros(numGenesTreat,1);
        for j = 1:numGenesTreat
            counts{k}(j) = sum(weights(strcmp(geneLists{k}{j},geneList)));
        end
    else
        % Just count how many times it appears in the agglomerated list
        % (i.e., how many drugs a gene appears in)
        counts{k} = [frequencyTable{:,2}];
    end
    
    % Sort by counts:
    [~,ix] = sort(counts{k},'descend');
    geneLists{k} = geneLists{k}(ix);
    counts{k} = counts{k}(ix);
    
    % Talk to me:
    fprintf(1,'%u genes are targeted by existing drugs for %s:\n',...
        numGenesTreat,whatDisease);
    for i = 1:numGenesTreat
        %fprintf(1,'%s, ',geneLists{k}{i});
    end
    fprintf(1,'\n');
end

%-------------------------------------------------------------------------------
% Process into a unique set of genes:
allGenes = unique(vertcat(geneLists{:}));
numGenes = length(allGenes);
fprintf(1,'%u genes assigned to at least one drug across our %u disorders\n',...
    numGenes,numDiseases);

%-------------------------------------------------------------------------------
% Construct a gene x disease table
indicatorMatrix = zeros(numGenes,numDiseases);
for i = 1:numGenes
    for j = 1:numDiseases
        theGeneIsHere = strcmp(allGenes{i},geneLists{j});
        if any(theGeneIsHere)
            indicatorMatrix(i,j) = counts{j}(theGeneIsHere);
        end
    end
end

%-------------------------------------------------------------------------------
% Convert to proportions:
propMatrix = indicatorMatrix;
for j = 1:numDiseases
    propMatrix(:,j) = propMatrix(:,j)/sum(propMatrix(:,j))*100;
end
meanRow = mean(propMatrix,2);

% Normalize the indicator matrix:
if normalizeWithinDrugs
    % Normalize so each disease gets the same total weight to distribute:
    for k = 1:numDiseases
        indicatorMatrix(:,k) = indicatorMatrix(:,k)/sum(indicatorMatrix(:,k));
    end
    fprintf(1,'Normalizing indicator matrix\n');
end

% Make a table:
indicatorTable = array2table(indicatorMatrix,'RowNames',allGenes,...
    'VariableNames',whatDiseases);
percIndicatorTable = array2table(propMatrix,'RowNames',allGenes,...
    'VariableNames',whatDiseases);

% Sort the table:
[~,ix] = sort(meanRow,'descend');
indicatorTable = indicatorTable(ix,:);
percIndicatorTable = percIndicatorTable(ix,:);

end
