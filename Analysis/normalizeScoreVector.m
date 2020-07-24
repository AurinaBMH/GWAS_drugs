function geneWeightsNorm = normalizeScoreVector(geneWeights, whatNorm)
if nargin < 2
    whatNorm = 1; 
end

r = ~isnan(geneWeights);
geneWeightsNorm = geneWeights;
geneWeightsNorm(r) = geneWeights(r)/norm(geneWeights(r),whatNorm);

end