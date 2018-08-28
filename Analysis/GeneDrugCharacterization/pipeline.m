function geneScores = pipeline(whatDisease)
% Pipeline for producing table characterizing individual genes
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Parse inputs:
%-------------------------------------------------------------------------------
if nargin < 1
    % GWAS hits from which disease?: 'all','SZP','ASD','ADHD','BIP','MDD'
    whatDisease = 'MDD';
end

%-------------------------------------------------------------------------------
% Load in default parameters:
%-------------------------------------------------------------------------------
params = SetDefaultParams();
doWeighted = params.doWeighted;
LDthreshold = params.LDthreshold;

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

%-------------------------------------------------------------------------------
% Import SNP, gene, disease, GWAS, LD annotations for a given GWAS study:
%-------------------------------------------------------------------------------
SNPAnnotationTable = SNPAnnotationImport(whatDisease,LDthreshold);

% Combine all mapped genes:
allMappedDiseaseGenes = unique(vertcat(SNPAnnotationTable.mappedGenes{:}));

% Add LD-genes:
onlyLDDiseaseGenes = unique(vertcat(SNPAnnotationTable.LDgenes{:}));
allLDDiseaseGenes = union(allMappedDiseaseGenes,onlyLDDiseaseGenes);

%===============================================================================
%======================CHARACTERIZATION MODULES=================================
%===============================================================================
% Now get different characterizations of each gene of interest to closeness
% to different data types.
%-------------------------------------------------------------------------------
geneScores = struct();
geneScores.params = params;

%-------------------------------------------------------------------------------
% SNP->gene distance on DNA (mapped SNP, or LD SNP to a gene)
%-------------------------------------------------------------------------------
geneScores.DNA = TellMeDNADistance(allUniqueGenes,SNPAnnotationTable);

%-------------------------------------------------------------------------------
% eQTL Characterization:
%-------------------------------------------------------------------------------
% eQTL_info = TellMeQTL(allDiseaseGenesMapped,allUniqueGenes);

%-------------------------------------------------------------------------------
% PPIN characterization:
%-------------------------------------------------------------------------------

% (i) binarized at zero evidence threshold:
numSteps = 3;
geneScores.PPI_mapped_th0 = TellMePPIInfo(allMappedDiseaseGenes,allUniqueGenes,false,0);
geneScores.PPI_LD_th0 = TellMePPIInfo(allLDDiseaseGenes,allUniqueGenes,false,0);

% (ii) binarized at an evidence threshold of 0.4:
numSteps = 4;
geneScores.PPI_mapped_th400 = TellMePPIInfo(allMappedDiseaseGenes,allUniqueGenes,false,400,numSteps);
geneScores.PPI_LD_th400 = TellMePPIInfo(allLDDiseaseGenes,allUniqueGenes,false,400,numSteps);

% (iii) weighted:
% numSteps = 4;
% % Just mapped disease genes (genes with a GWAS SNP in them):
% geneScores.PPI_mapped_weighted = TellMePPIInfo(allMappedDiseaseGenes,allUniqueGenes,true,numSteps);
% % Include genes LD to GWAS SNPs:
% geneScores.PPI_LD_weighted = TellMePPIInfo(allLDDiseaseGenes,allUniqueGenes,true,numSteps);

%-------------------------------------------------------------------------------
% AHBA gene coexpression:
%-------------------------------------------------------------------------------
% For mapped SNPs:
geneScores.AllenMeanCoexpMapped = TellMeAllenCoexp(allUniqueGenes,allMappedDiseaseGenes);
% Including LD SNPs:
geneScores.AllenMeanCoexpLD = TellMeAllenCoexp(allUniqueGenes,allLDDiseaseGenes);

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
