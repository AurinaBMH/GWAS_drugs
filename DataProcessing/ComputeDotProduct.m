function rho = ComputeDotProduct(v1,v2,doShuffle,minGoodProp,whatNorm)
% Compute the weighted sum of two vectors, v1 & v2
if nargin < 3
    doShuffle = false;
end
if nargin < 4
    % Require fewer than 40% bad values after matching indices
    minGoodProp = 0.4;
end
if nargin < 5
    whatNorm = 2;
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Checks and processing
%-------------------------------------------------------------------------------
% Check if too many bad values:
r = ~isnan(v1) & ~isnan(v2);
if mean(r) < minGoodProp
    rho = NaN;
    return
end

if any(v1 < 0) || any(v2 < 0)
    warning('Weight vectors contain negative values')
end

% Shuffle:
if doShuffle
    % Shuffle weights on v1:
    v1 = v1(randperm(length(v1)));
end

% Filter to fixed sparsity, on top X scores:
% X = 200;
% [~,ix] = sort(v1,'descend');
% v1(v1<v1(ix(X))) = 0;
% [~,ix] = sort(v2,'descend');
% v2(v2<v2(ix(X))) = 0;

%-------------------------------------------------------------------------------
% Compute dot product, ignoring dimensions containing NaNs
%-------------------------------------------------------------------------------
v1(isnan(v1)) = 0;
v2(isnan(v2)) = 0;
v1 = v1/norm(v1,whatNorm);
v2 = v2/norm(v2,whatNorm);
rho = sum(v1.*v2);

% v1_good = v1(r);
% v2_good = v2(r);
% if doShuffle
%     % Shuffle weights on v1:
%     v1_good = v1_good(randperm(length(v1_good)));
% end
%
% % Normalize v1 and v2 to unit vectors:
% if all(v1_good > 0)
%     v1_good = v1_good/norm(v1_good,2);
% end
% if all(v2_good > 0)
%     v2_good = v2_good/norm(v2_good,2);
% end

% Compute the score:
% rho = sum(v1_good.*v2_good);
end
