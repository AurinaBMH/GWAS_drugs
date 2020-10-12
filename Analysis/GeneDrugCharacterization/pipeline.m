function geneScores = pipeline(DISORDERlist, whatDisease, whatThreshold)

% use all mapping methods in the file
mappingMethods = fieldnames(DISORDERlist); 

% Pipeline for producing table characterizing individual genes
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Parse inputs:
%-------------------------------------------------------------------------------
if nargin < 1
    % GWAS hits from which disease?: 'all','SZP','ASD','ADHD','BIP','MDD'
    whatDisease = 'MDD2';
end

%-------------------------------------------------------------------------------
% Load in default parameters:
%-------------------------------------------------------------------------------
params = SetDefaultParams();
doWeighted = params.doWeighted;
geneScore = params.geneScore; 
whatDrugTargets = params.whatDrugTargets; 

%===============================================================================
% LOAD DATA
%===============================================================================

%-------------------------------------------------------------------------------
% Restrict our characterization to genes associated with (any) drug-treatment:
%-------------------------------------------------------------------------------
% List of all genes with drug targets in Drugbank are in **gene_ATC_matrix.csv**
indicatorTable = ImportTreatmentLists(false, whatDrugTargets);
allUniqueGenes = indicatorTable.Properties.RowNames;
numUniqueGenes = length(allUniqueGenes);
fprintf(1,'Analyzing %u genes that have drug targets in our list\n',numUniqueGenes);

%-------------------------------------------------------------------------------
% Get genes for a given GWAS study for different mapping methods:
%-------------------------------------------------------------------------------

% Directly mapped gened from standard MAGMA analysis will be used as inputs to PPI network and gene coexpression calculations
% Set threshold to include diretcly mapped genes to PPI network
% According to Zac's recommendation this could be a 0.05/number of genes identified through the MAGMA analysis; 
listGENESmapped = DISORDERlist.MAGMAdefault.(whatDisease);
listGENESeQTLbrain = DISORDERlist.eQTLbrain.(whatDisease);

switch whatThreshold
    case 'BF'
        pThr_m = 0.05/size(listGENESmapped,1); % Bonf correction for the number of genes in the list
        pThr_e = 0.05/size(listGENESeQTLbrain,1); % Bonf correction for the number of genes in the list
        allMappedDiseaseGenes = listGENESmapped.GENENAME(listGENESmapped.P<pThr_m);
        alleQTLbrainDiseaseGenes = listGENESeQTLbrain.GENENAME(listGENESeQTLbrain.P<pThr_e);
        
    case 'FDR'
        pFDR_m = mafdr(listGENESmapped.P, 'BHFDR', true);
        pFDR_e = mafdr(listGENESeQTLbrain.P, 'BHFDR', true);
        allMappedDiseaseGenes = listGENESmapped.GENENAME(pFDR_m<0.05);
        alleQTLbrainDiseaseGenes = listGENESeQTLbrain.GENENAME(pFDR_e<0.05);
end


%===============================================================================
%======================CHARACTERIZATION MODULES=================================
%===============================================================================
% Now get different characterizations of each gene of interest to closeness
% to different data types.
%-------------------------------------------------------------------------------
geneScores = struct();
geneScores.params = params;

geneWeights = {'ZSTAT', 'P', 'NSNPS', 'NSNPSnorm'}; 
%-------------------------------------------------------------------------------
% SNP->gene distance on DNA (mapped SNP, or LD SNP to a gene)
%-------------------------------------------------------------------------------
for m=1:length(mappingMethods)
    
    for s=1:length(geneWeights)
        % give an empty vector
        scoreVector = zeros(length(allUniqueGenes),1);
        
        mapping = mappingMethods{m};
        mapTABLE = DISORDERlist.(mapping).(whatDisease);
        % select genes that are drug targets
        [~, INDput, INDtake] = intersect(allUniqueGenes, mapTABLE.GENENAME);
        
        % give the selected score to all mentioned genes and keep 0 to all others
        switch geneWeights{s}
            case 'NSNPSnorm'
                % give a score of NSNPS normalised to gene length
                geneLength = mapTABLE.STOP-mapTABLE.START;
                scoreVector(INDput) = mapTABLE.NSNPS(INDtake)./geneLength(INDtake);
            case 'P'
                % take the log10p to indicate p-values - higher values are "better"
                scoreVector(INDput) = -log10(mapTABLE.(geneWeights{s})(INDtake));
            otherwise
                scoreVector(INDput) = mapTABLE.(geneWeights{s})(INDtake);
        end
        
        geneScores.(mapping).(geneWeights{s}) = scoreVector;
        
    end
end

%-------------------------------------------------------------------------------
% PPIN characterization:
%-------------------------------------------------------------------------------
% Considering:
% (i) Just mapped disease genes (genes with a GWAS SNP in them) 
% - these will be from MAGMA analysis
% (ii) Include genes that are brain eQTLs 

% (*) binarized at zero evidence threshold:
numSteps = 3;
geneScores.PPI_mapped_th0 = TellMePPIInfo(allMappedDiseaseGenes,allUniqueGenes,false,0,numSteps);
geneScores.PPI_eQTLbrain_th0 = TellMePPIInfo(alleQTLbrainDiseaseGenes,allUniqueGenes,false,0,numSteps);

% (*) binarized at an evidence threshold of 0.4:
numSteps = 3;
geneScores.PPI_mapped_th400 = TellMePPIInfo(allMappedDiseaseGenes,allUniqueGenes,false,400,numSteps);
geneScores.PPI_eQTLbrain_th400 = TellMePPIInfo(alleQTLbrainDiseaseGenes,allUniqueGenes,false,400,numSteps);

% (*) binarized at an evidence threshold of 0.6:
numSteps = 3;
geneScores.PPI_mapped_th600 = TellMePPIInfo(allMappedDiseaseGenes,allUniqueGenes,false,600,numSteps);
geneScores.PPI_eQTLbrain_th600 = TellMePPIInfo(alleQTLbrainDiseaseGenes,allUniqueGenes,false,600,numSteps);

% (*) Binarized at an evidence threshold of 0.9:
numSteps = 3;
geneScores.PPI_mapped_th900 = TellMePPIInfo(allMappedDiseaseGenes,allUniqueGenes,false,900,numSteps);
geneScores.PPI_eQTLbrain_th900 = TellMePPIInfo(alleQTLbrainDiseaseGenes,allUniqueGenes,false,900,numSteps);

%-------------------------------------------------------------------------------
% AHBA gene coexpression:
%-------------------------------------------------------------------------------
% For mapped genes: 
geneScores.AllenMeanCoexpMapped = TellMeAllenCoexp(allUniqueGenes,allMappedDiseaseGenes);
% Including barin eQTL genes:
geneScores.AllenMeanCoexpeQTLbrain = TellMeAllenCoexp(allUniqueGenes,alleQTLbrainDiseaseGenes);

%===============================================================================
% Assimilate results
%===============================================================================

% Now make a table
geneScores.gene = allUniqueGenes;

end
