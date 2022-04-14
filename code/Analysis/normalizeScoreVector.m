function geneWeightsNorm = normalizeScoreVector(geneWeights, whatNorm)
if nargin < 2
     whatNorm = 2; % as was in the original script 
%     % 1-sum of the absolute values of the vector elements.
%     % 2-vector magnitude or Euclidean length of the vector.
end

r = ~isnan(geneWeights);
geneWeightsNorm = geneWeights;
geneWeightsNorm(r) = geneWeights(r)/norm(geneWeights(r),whatNorm);

end