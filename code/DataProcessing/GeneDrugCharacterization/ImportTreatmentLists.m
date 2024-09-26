function [indicatorTable,percIndicatorTable, dataTable, allDrugs] = ImportTreatmentLists(normalizeWithinDrugs, whatDrugTargets, whatTargets)
% Import information on gene targets for psychiatric conditions
%-------------------------------------------------------------------------------
% Input parameters:
%-------------------------------------------------------------------------------
if nargin < 1
    normalizeWithinDrugs = true;
end
if nargin < 2
    whatDrugTargets = '2024'; 
    params = SetDefaultParams();
    whatTargets = params.whatTargets; 
    % 2020 - uses automated AA version
end
if nargin < 3
    params = SetDefaultParams();
    whatTargets = params.whatTargets; 
end


if normalizeWithinDrugs
    fprintf(1,'Weighting genes equally within a drug...\n');
end

%-------------------------------------------------------------------------------
% Whether to treat each drug as equally important, and each gene targeted by a
% drug as equally important for the efficacy of that drug:

params = SetDefaultParams();

%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------

switch whatDrugTargets
    case '2020'
        % use drug targets assigned automatically by AA in 08/2020
        % it takes ~10s to run, so load the pre-computed data here
        % dataTable = give_drugTargets();
        
        fileName = sprintf('DataOutput_2022/drugTargets_2020_%s_drugbank.mat', whatTargets);
        load(fileName, 'dataTable');
        
        whatDiseases = params.whatDiseases_Treatment_ALL; 
        numDiseases = length(whatDiseases);

        
    case '2024'
        
        fileName = sprintf('DataOutput_2024/drugTargets_2024_%s_drugbank.mat', whatTargets);
        load(fileName, 'dataTable');
        
        whatDiseases = params.whatDiseases_Treatment_ALL; 
        numDiseases = length(whatDiseases);
   
end

% get all drugs with active gene targets from DrugBank
allDrugs = get_allDrugBank_targets(whatTargets); 
% make a table of all mentioned drugs with their targets keeping only unique ones 
% DT = cell(numDiseases,1); 
% for kk=1:numDiseases
%    DT{kk} = dataTable.(whatDiseases{kk}); 
% end 
% 
% allDrugs = vertcat(DT{:}); 
% % remove drugs with non-existent targets
% allDrugs(strcmp(string(allDrugs.Target), ""),:) = [];  
% % select only unique drug names
% [~, INDunique] = unique(allDrugs.Name); 
% allDrugs = allDrugs(INDunique,:); 

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
% unique set of genes now consists of all available genes as active DrugBank targets; 
% allGenes = unique(vertcat(geneLists{:}));
allGenes = strjoin(allDrugs.Target,', ');
allGenes = unique(split(allGenes, ', '));
allGenes = allGenes(~cellfun('isempty',allGenes));  
numGenes = length(allGenes);
fprintf(1,'There are %u drug target genes in the DrugBank\n', numGenes);

%-------------------------------------------------------------------------------
% Construct a gene x disease table
indicatorMatrix = zeros(numGenes,numDiseases);
for i = 1:numGenes
    for j = 1:numDiseases
        theGeneIsHere = strcmpi(allGenes{i},geneLists{j});
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
% [~,ix] = sort(meanRow,'descend');
% indicatorTable = indicatorTable(ix,:);
% percIndicatorTable = percIndicatorTable(ix,:);

end
