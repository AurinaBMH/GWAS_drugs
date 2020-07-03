% binary full
[AdjPPI,geneNames] = PPINImport(false,0,'HGNCmatch'); 

% binary >400 evidence score
[AdjPPI,geneNames] = PPINImport(false,400,'HGNCmatch');

% weighted
[AdjPPI,geneNames] = PPINImport(true);

% compute distances between genes 
distMatrix = ComputePPIDist(400,false);