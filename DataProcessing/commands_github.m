% binary full
[AdjPPI,geneNames] = PPINImport(false,0,'HGNCmatch'); 

% binary >400 evidence score
[AdjPPI,geneNames] = PPINImport(false,400,'HGNCmatch');

% weighted
[AdjPPI,geneNames] = PPINImport(true);

% compute distances between genes 
distMatrix = ComputePPIDist(0,false);
distMatrix = ComputePPIDist(400,false);
distMatrix = ComputePPIDist(600,false);
distMatrix = ComputePPIDist(900,false);

distMatrix = ComputePPIDist([],true); - % this takes ages, didn't finish in 3 days. 

% now result tables are generated for each disorder
GenerateResultsTables; 

% loop over all mapping methods and save the results to file, so can compare; 
% load ADHD as example to select mapping methods
load('resultsTable_ADHD_FDR.mat', 'geneScores')
similarityTypes = setdiff(fieldnames(geneScores), {'gene', 'params'}); 
whatThreshold = 'BF'; 

for t=1:length(similarityTypes)
    
    if contains(similarityTypes{t},'PPI')
        whatProperty = 'percPPIneighbors1';

    elseif contains(similarityTypes{t},'Allen')
        whatProperty = 'r';
    else
        whatProperty = 'ZSTAT';
    end
    
    DistinguishingCharBar(similarityTypes{t},whatProperty, whatThreshold)
%     figureName = sprintf('figures/GWASdrug_%s_%s_%s', similarityTypes{t},whatProperty, whatThreshold);
%     print(gcf,figureName,'-dpng','-r300');

end



