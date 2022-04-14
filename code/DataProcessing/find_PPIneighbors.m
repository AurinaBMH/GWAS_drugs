function PPI_neighborsK = find_PPIneighbors(PPIN, geneList, numSteps)

numGenesChar = length(geneList);
PPI_neighbors_gene = cell(numSteps, numGenesChar); 

for i = 1:numGenesChar
    
    gene_i = geneList{i};
    
    % First try matching the gene:
    PPI_index = find(strcmp(PPIN.geneNames,gene_i));
    
    if isempty(PPI_index)
        fprintf('%s is not represented in the PPI network\n',gene_i)
        continue;
    end
    
    %-------------------------------------------------------------------------------
    % Get the k-step neighbors:
    for k = 1:numSteps
        % (i) Don't distinguish neighbors by pathlength (take union):
        iskStepNeighbor = (PPIN.distMatrix(PPI_index,:) <= k);
        
        % Sum of the number of labeled genes that are within k steps (maximum score is the number of labeled genes)
        numNeighbors = sum(iskStepNeighbor);
        if numNeighbors==0
            PPI_neighbors_gene{k,i} = NaN;
        else
            PPI_neighbors_gene{k,i} = PPIN.geneNames(iskStepNeighbor);
            
        end
    end

end

PPI_neighborsK = cell(numSteps); 

for k = 1:numSteps
    PPI_neighbors = PPI_neighbors_gene(k,:);
    PPI_neighborsK{k} = unique(cat(1, PPI_neighbors{:}));
end

% concatonate a list of k-step neighbors to get a single list of neighbors
end