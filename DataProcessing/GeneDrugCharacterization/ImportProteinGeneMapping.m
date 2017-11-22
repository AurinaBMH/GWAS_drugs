function dataTable = ImportProteinGeneMapping(allUniqueGenes,PPINGeneNames)
%-------------------------------------------------------------------------------
% Idea is to match HGNC gene names to names given in the PPIN network
%-------------------------------------------------------------------------------

% Preliminaries:
numGenes = length(allUniqueGenes);
fprintf(1,'%u genes to match\n',numGenes);
onlyShowSuccess = true; % only give feedback when works
% keyboard

%-------------------------------------------------------------------------------
%% Import data from text file.
%-------------------------------------------------------------------------------
%% Initialize variables.
filename = 'PPI_forIDmatching_final.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
fileID = fopen(filename,'r');
%% Read columns of data according to the format:
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
%% Create output variable
dataTable = table(dataArray{1:end-1}, 'VariableNames', {'STRINGids','HGNC_1','HGNC_2','HGNC_3','HGNC_4','HGNC_5','HGNC_6','HGNC_7','HGNC_8','HGNC_9','HGNC_10','HGNC_11','HGNC_12','HGNC_13','HGNC_14','HGNC_15','HGNC_16'});
%-------------------------------------------------------------------------------

% Now we can try to match:
allUniqueProteins = cell(numGenes,1);
for i = 1:numGenes
    geneName = allUniqueGenes{i};
    % Try to match this gene
    matchMe = strcmp(geneName,PPINGeneNames);
    numMatches = sum(matchMe);
    if numMatches==0
        if ~onlyShowSuccess
            fprintf(1,'No direct PPIN match for %s\n',geneName);
        end
        % Look for matches in the dataTable
        matchWhat = 0;
        for j = 1:16
            colName = sprintf('HGNC_%u',j);
            matchMeHGNC = strcmp(geneName,dataTable.(colName));
            numMatchesHGNC = sum(matchMeHGNC);
            if numMatchesHGNC > 0
                matchWhat = find(matchMeHGNC);
                % fprintf(1,'No direct match for %s but matched to %s!\n',geneName,dataTable.(colName){matchWhat});
                break
            end
        end
        if matchWhat > 0
            % Found a match
            % see if any match the PPINGeneNames
            theSynonyms = GiveMeSynonyms(dataTable,matchWhat);
            matchPPIN = ismember(theSynonyms,PPINGeneNames);
            if any(matchPPIN)
                fprintf(1,'~~~OMG it actually happened! We got a match after looking at synonyms!\n');
                fprintf(1,'Using %s (matches) instead of %s (doesn''t match)\n',...
                                theSynonyms{matchPPIN},geneName);
                allUniqueProteins{i} = theSynonyms{matchPPIN};
            else
                if ~onlyShowSuccess
                    fprintf(1,'No HGNC synonyms match :/\n');
                end
            end
        else
            if ~onlyShowSuccess
                fprintf(1,'%s not found in HGNC ID-matching table...??\n',geneName);
            end
        end

    elseif numMatches==1
        % fprintf(1,'Too easy boyz -- %s matches already aight\n',geneName);
        allUniqueProteins{i} = geneName;
    end
end

%-------------------------------------------------------------------------------
function theSynonyms = GiveMeSynonyms(dataTable,indx)
    % Return a string at indx
    theSynonyms = cell(1,1);
    for j = 1:16
        colName = sprintf('HGNC_%u',j);
        tableValue = dataTable.(colName){indx};
        if strcmp(tableValue,'NA')
            break
        else
            theSynonyms{j} = tableValue;
        end
    end
end

end
