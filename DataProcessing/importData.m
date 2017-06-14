eQTLedges = importEdges();
eQTLidentifier = importIdentifier();
eQTLproteinnames = importProteinNames();

numProteins = length(eQTLproteinnames);

fprintf(1,'Constructing a %ux%u protein-protein interaction network\n',...
                    numProteins,numProteins);

Adj = zeros(numProteins,numProteins);
