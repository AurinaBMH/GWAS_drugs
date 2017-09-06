function importData(doEQTL)

if nargin < 1
    doEQTL = true;
end

%-------------------------------------------------------------------------------
% Do the importing from csv files:
%-------------------------------------------------------------------------------
fprintf(1,'Importing data from csv files\n');
if doEQTL
    geneIdentifier = importIdentifier_eQTL();
    proteinNames = importProteinNames();
    PPI_edges = importEdges_eQTL();
else
    geneIdentifier = ImportIdentifierMapped();
    proteinNames = geneIdentifier.Name;
    PPI_edges = importEdgesMapped();
end
numProteins = length(proteinNames);

%-------------------------------------------------------------------------------
% Interaction matrix
%-------------------------------------------------------------------------------
numInteractions = height(PPI_edges);

fprintf(1,['Constructing a %ux%u protein-protein interaction network',...
                ' containing %u interactions...\n'],...
                    numProteins,numProteins,numInteractions);

Adj = zeros(numProteins,numProteins);

for k = 1:numInteractions
    ii = strcmp(proteinNames,PPI_edges.InteractorA(k));
    jj = strcmp(proteinNames,PPI_edges.InteractorB(k));
    Adj(ii,jj) = 1;
end
Adj = (Adj | Adj');

%-------------------------------------------------------------------------------
% Process node annotations
%-------------------------------------------------------------------------------
% Agglomerate data in geneIdentifier
% proteinMetaData = table();
isGWAS = ismember(proteinNames,unique(geneIdentifier.Name(geneIdentifier.GWAS)));
isLD = ismember(proteinNames,unique(geneIdentifier.Name(geneIdentifier.LD)));
if doEQTL
    isNeighbor = ismember(proteinNames,unique(geneIdentifier.Name(geneIdentifier.Partners)));
else
    % No neighbors in the matched version:
    isNeighbor = false(numProteins,1);
end
fprintf(1,'%u GWAS, %u LD, %u neighbor\n',sum(isGWAS),sum(isLD),sum(isNeighbor));


%-------------------------------------------------------------------------------
% Save to file
%-------------------------------------------------------------------------------

if doEQTL
    fileName = 'processedData_eQTL.mat';
else
    fileName = 'processedData_Matched.mat';
end
save(fullfile('Data',fileName),'Adj','proteinNames');
fprintf(1,'Saved processed data to ''%s''\n',fileName);
