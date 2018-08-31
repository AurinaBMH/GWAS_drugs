function D = distance_bin(A)
%DISTANCE_BIN       Distance matrix
%
%   D = distance_bin(A);
%
%   The distance matrix contains lengths of shortest paths between all
%   pairs of nodes. An entry (u,v) represents the length of shortest path
%   from node u to node v. The average shortest path length is the
%   characteristic path length of the network.
%
%   Input:      A,      binary directed/undirected connection matrix
%
%   Output:     D,      distance matrix
%
%   Notes:
%       Lengths between disconnected nodes are set to Inf.
%       Lengths on the main diagonal are set to 0.
%
%   Algorithm: Algebraic shortest paths.
%
%   Mika Rubinov, U Cambridge
%   Jonathan Clayden, UCL
%   2007-2013

% Modification history:
% 2007: Original (MR)
% 2013: Bug fix, enforce zero distance for self-connections (JC)

%-------------------------------------------------------------------------------
% Checks on input matrix:
if any(A(:) < 0)
    error('Adjacency matrix contains negative values')
end
if issparse(A)
    fprintf(1,'Converting sparse to a full matrix :O\n');
    A = full(A);
end
A = double(A > 0); % binarize and convert to double format
if any(diag(A) > 0); % Self connections lead to a never-ending while loop
    warning('%u self-connections?? Removed.',sum(diag(A) > 0))
    A(logical(eye(size(D)))) = 0;
end
%-------------------------------------------------------------------------------

l = 1;             % path length
Lpath = A;         % matrix of paths l
D = A;             % distance matrix

% Iterate through path lengths:
Idx = true;
while any(Idx(:))
    fprintf(1,'Path length %u...\n',l);
    l = l + 1;
    Lpath = Lpath*A;
    Idx = (Lpath~=0) & (D==0);
    fprintf(1,'%u pairs of nodes have a path length of %u...\n',sum(Idx(:)),l);
    D(Idx) = l;
end

% Clean up:
D(~D) = Inf;                      % assign inf to disconnected nodes
D(logical(eye(size(D)))) = 0;     % clear diagonal

end
