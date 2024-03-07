% check how many genes survive thresholding

load('GWAS_disordersMAGMA_2024.mat')

mappings = {'MAGMAdefault','eQTLbrain'}; 
disorders = fieldnames(DISORDERlist.MAGMAdefault); 
        
for i=1:length(mappings)
    for j=1:length(disorders)
        
        listGENESmapped = DISORDERlist.(mappings{i}).(disorders{j}); 
        pThr_m = 0.05/size(listGENESmapped,1); % Bonf correction for the number of genes in the list
        
        allMappedDiseaseGenes = listGENESmapped.GENENAME(listGENESmapped.P<pThr_m);

        fprintf('%s %s: %d genes\n', mappings{i}, disorders{j}, length(allMappedDiseaseGenes));
        
    end
end



% visualise PPI network distances at different thresholds
thrs = [400,600,900]; 
figure;
for k=1:length(thrs)
    
load(sprintf('PPI_HGNC_Dist_th%d.mat', thrs(k)))
subplot(1,3,k); 
numDis = length(find(isinf(distMatrix(:,1))));
numGenes = size(distMatrix,1); 
histogram(maskuHalf(distMatrix)); title(sprintf('PPI distances thr%d, %d/%d genes disconnected', thrs(k), numDis, numGenes))

end


% check the number of neighbors in PPI network for one disorder
load('resultsTable_BIP2.mat')

% Count disease genes that are X-step neighbors on the PPI network
figure; 
subplot(1,3,1); histogram(geneScores.PPI_eQTLbrain_th600.numPPIneighbors1); title('numPPIneighbors1')
subplot(1,3,2); histogram(geneScores.PPI_eQTLbrain_th600.numPPIneighbors2); title('numPPIneighbors2')% this has right-skewed distribution
subplot(1,3,3); histogram(geneScores.PPI_eQTLbrain_th600.numPPIneighbors3); title('numPPIneighbors3')



figure; histogram(geneScores.PPI_eQTLbrain_th0.numPPIneighbors1); 
figure; histogram(geneScores.PPI_eQTLbrain_th0.numPPIneighbors2);
figure; histogram(geneScores.PPI_eQTLbrain_th0.numPPIneighbors3); 
        

figure; histogram(geneScores.PPI_eQTLbrain_th400.numPPIneighbors1); 
figure; histogram(geneScores.PPI_eQTLbrain_th400.numPPIneighbors2); % this has normal distribution
figure; histogram(geneScores.PPI_eQTLbrain_th400.numPPIneighbors3); 



figure; histogram(geneScores.PPI_eQTLbrain_th900.numPPIneighbors1); 
figure; histogram(geneScores.PPI_eQTLbrain_th900.numPPIneighbors2); % this has right-skewed distribution
figure; histogram(geneScores.PPI_eQTLbrain_th900.numPPIneighbors3); 
figure; histogram(geneScores.PPI_eQTLbrain_th900.numPPIneighbors4); 

