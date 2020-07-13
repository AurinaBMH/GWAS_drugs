function geneScores = pipeline(DISORDERlist, whatDisease)

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

%===============================================================================
% LOAD DATA
%===============================================================================

%-------------------------------------------------------------------------------
% Restrict our characterization to genes associated with (any) drug-treatment:
%-------------------------------------------------------------------------------
% List of all genes with drug targets in Drugbank are in **gene_ATC_matrix.csv**
indicatorTable = ImportTreatmentLists(false);
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

pThr_m = 0.05/size(listGENESmapped,1); % Bonf correction for the number of genes in the list 
pThr_e = 0.05/size(listGENESeQTLbrain,1); % Bonf correction for the number of genes in the list 

allMappedDiseaseGenes = listGENESmapped.GENENAME(listGENESmapped.P<=pThr_m); 
alleQTLbrainDiseaseGenes = listGENESeQTLbrain.GENENAME(listGENESeQTLbrain.P<=pThr_e); 

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
numSteps = 5;
geneScores.PPI_mapped_th0 = TellMePPIInfo(allMappedDiseaseGenes,allUniqueGenes,false,0,numSteps);
geneScores.PPI_eQTLbrain_th0 = TellMePPIInfo(alleQTLbrainDiseaseGenes,allUniqueGenes,false,0,numSteps);

% (*) binarized at an evidence threshold of 0.4:
numSteps = 6;
geneScores.PPI_mapped_th400 = TellMePPIInfo(allMappedDiseaseGenes,allUniqueGenes,false,400,numSteps);
geneScores.PPI_eQTLbrain_th400 = TellMePPIInfo(alleQTLbrainDiseaseGenes,allUniqueGenes,false,400,numSteps);

% (*) binarized at an evidence threshold of 0.6:
numSteps = 6;
geneScores.PPI_mapped_th600 = TellMePPIInfo(allMappedDiseaseGenes,allUniqueGenes,false,600,numSteps);
geneScores.PPI_eQTLbrain_th600 = TellMePPIInfo(alleQTLbrainDiseaseGenes,allUniqueGenes,false,600,numSteps);

% (*) Binarized at an evidence threshold of 0.9:
numSteps = 6;
geneScores.PPI_mapped_th900 = TellMePPIInfo(allMappedDiseaseGenes,allUniqueGenes,false,900,numSteps);
geneScores.PPI_eQTLbrain_th900 = TellMePPIInfo(alleQTLbrainDiseaseGenes,allUniqueGenes,false,900,numSteps);

% (*) weighted: - weighted PPI distances are not calculated for now
% numSteps = 6;

%geneScores.PPI_mapped_weighted = TellMePPIInfo(allMappedDiseaseGenes,allUniqueGenes,true,numSteps);
%geneScores.PPI_eQTLbrain_weighted = TellMePPIInfo(alleQTLbrainDiseaseGenes,allUniqueGenes,true,numSteps);

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
gene = allUniqueGenes;
% resultsTable = table(gene,geneScores.PPI_mapped_weighted,numGWASMapped,numLDSNPs,percPPIneigh1Mapped,...
%                         percPPIneigh1LD,AllenMeanCoexpMapped,AllenMeanCoexpLD);
% % ,numSNPGenes,numEGenes_LD,numSNPGenes_LD,numLDeGeneseQTL,numLDSNPGeneseQTL); % matchingDrugsString
%
% % Sort by column, then by column, etc. in ordered hierarchy:
% resultsTable = sortrows(resultsTable,{'numGWASMapped','numLDSNPs','percPPIneigh1Mapped',...
%                 'percPPIneigh1LD','AllenMeanCoexpMapped',...
%                 'AllenMeanCoexpLD'},'descend','MissingPlacement','last');


%-------------------------------------------------------------------------------
% Display just with custom columns
% customColumns = {'gene','numGWASMapped','numLDSNPs','percPPIneigh1Mapped',...
%                 'percPPIneigh1LD','AllenMeanCoexpMapped'}; % ,'matchingDrugsString'
% display(resultsTable(1:40,ismember(resultsTable.Properties.VariableNames,customColumns)));

%===============================================================================
% for i = 1:numUniqueGenes
%     gene_i = allUniqueGenes{i};
%     fprintf(1,'[%u/%u]: %s (%s)\n',i,numUniqueGenes,gene_i,protein_i);
%
%     %-------------------------------------------------------------------------------
%     % Preliminaries:
%     %-------------------------------------------------------------------------------
%     % (I assume this can be comprehensive given the data provided??? Maybe not??)
%     theLDgenes = GiveMeLDGenes(gene_i,SNPGeneMap,LDRelateTable,allDiseaseSNPs);
%     fprintf(1,'%u genes LD to the target\n',length(theLDgenes));
%
%
%     %-------------------------------------------------------------------------------
%     % Drugs, classes
%     %-------------------------------------------------------------------------------
%     % Match to drugs using geneDrugTable
%     % assign class using drugClassTable
%     % matchingDrugs = geneDrugTable.drugName(strcmp(geneDrugTable.geneName,gene_i));
%     % assignClass = @(druggy)drugClassTable.whatClass(strcmp(drugClassTable.drugName,druggy));
%     % matchingDrugsClass = cellfun(@(x)sprintf('%s(%s)',x,char(assignClass(x))),...
%     %                     matchingDrugs,'UniformOutput',false);
%     % matchingDrugsString{i} = BF_cat(matchingDrugsClass);
% end


%-------------------------------------------------------------------------------
% Some basic stats as user info to screen:
% fprintf(1,'%u genes have mapped and LD\n',sum(resultsTable.numLDSNPs > 0 & resultsTable.numGWASMapped > 0));
% fprintf(1,'%u genes have mapped but no LD\n',sum(resultsTable.numLDSNPs==0 & resultsTable.numGWASMapped > 0));
% fprintf(1,'%u genes have LD but no mapped\n',sum(resultsTable.numLDSNPs > 0 & resultsTable.numGWASMapped==0));
% fprintf(1,'%u genes have eGene annotations\n',sum(resultsTable.numEGenes > 0));
% fprintf(1,'%u genes have SNPgene annotations\n',sum(resultsTable.numSNPGenes > 0));
% fprintf(1,'%u genes have eGene-LD annotations\n',sum(resultsTable.numEGenes_LD > 0));
% fprintf(1,'%u genes have SNPgene-LD annotations\n',sum(resultsTable.numSNPGenes_LD > 0));
% fprintf(1,'%u genes have LD-eGene annotations\n',sum(resultsTable.numEGenes_LD > 0));
% fprintf(1,'%u genes have LD-SNPgene annotations\n',sum(resultsTable.numSNPGenes_LD > 0));

%-------------------------------------------------------------------------------
% Display full table:
% display(resultsTable(1:60,:));



end
