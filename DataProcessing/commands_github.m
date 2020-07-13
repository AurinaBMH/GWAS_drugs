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
GenerateResultsTables


