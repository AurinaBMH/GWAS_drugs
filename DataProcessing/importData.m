%-------------------------------------------------------------------------------
% Do the importing from csv files:
%-------------------------------------------------------------------------------
fprintf(1,'Importing data from csv files\n');
eQTLedges = importEdges();
eQTLidentifier = importIdentifier();
eQTLproteinnames = importProteinNames();

%-------------------------------------------------------------------------------
% Interaction matrix
%-------------------------------------------------------------------------------
numInteractions = height(eQTLedges);
numProteins = length(eQTLproteinnames);

fprintf(1,['Constructing a %ux%u protein-protein interaction network',...
                ' containing %u interactions...\n'],...
                    numProteins,numProteins,numInteractions);

Adj = zeros(numProteins,numProteins);

for k = 1:numInteractions
    ii = strcmp(eQTLproteinnames,eQTLedges.InteractorA(k));
    jj = strcmp(eQTLproteinnames,eQTLedges.InteractorB(k));
    Adj(ii,jj) = 1;
end
Adj = (Adj | Adj');

%-------------------------------------------------------------------------------
% Process node annotations
%-------------------------------------------------------------------------------
% Agglomerate data in eQTLidentifier
% proteinMetaData = table();
isGWAS = ismember(eQTLproteinnames,unique(eQTLidentifier.Name(eQTLidentifier.GWAS)));
isLD = ismember(eQTLproteinnames,unique(eQTLidentifier.Name(eQTLidentifier.LD)));
isNeighbor = ismember(eQTLproteinnames,unique(eQTLidentifier.Name(eQTLidentifier.Partners)));
fprintf(1,'%u GWAS, %u LD, %u neighbor\n',sum(isGWAS),sum(isLD),sum(isNeighbor));


%-------------------------------------------------------------------------------
% Save to file
%-------------------------------------------------------------------------------
save(fullfile('Data','processedData.mat'),'Adj','eQTLproteinnames');
