function uniqueInd = countUnique(filter1,filter2,Adj)
% Idea is to count the number of unique genes (columns/rows) involved in a
% given type of interaction
%-------------------------------------------------------------------------------
N = size(Adj,1);

filterComb = repmat(filter1,1,N) & repmat(filter2,1,N)';
filterComb = (filterComb|filterComb');
filterConn = (Adj & filterComb);
uniqueInd = any(filterConn,2);

end
