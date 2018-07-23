function resultsTable = pipeline(whatDisease)
% Pipeline for producing table characterizing individual genes
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Parse inputs:
%-------------------------------------------------------------------------------
if nargin < 1
    % GWAS hits from which disease?: 'all','SZP','ASD','ADHD','BIP','MDD'
    whatDisease = 'all';
end

%-------------------------------------------------------------------------------
% Load in default parameters:
%-------------------------------------------------------------------------------
params = SetDefaultParams();
doWeighted = params.doWeighted;
PPINevidenceThreshold = params.PPINevidenceThreshold;
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
% to different data types:
[numGWASMapped,numLDSNPs] = TellMeDNADistance(allUniqueGenes,SNPAnnotationTable);

%-------------------------------------------------------------------------------
% eQTL Characterization:
%-------------------------------------------------------------------------------
% eQTL_info = TellMeQTL(allDiseaseGenesMapped,allUniqueGenes);

%-------------------------------------------------------------------------------
% PPIN characterization:
%-------------------------------------------------------------------------------
% Just mapped disease genes (genes with a GWAS SNP in them):
[numPPIneigh1Mapped,percPPIneigh1Mapped,meanPPIDistMapped] = TellMePPIInfo(PPINevidenceThreshold,allMappedDiseaseGenes,allUniqueGenes);

% Just genes LD to GWAS SNPs:
[numPPIneigh1LD,percPPIneigh1LD,meanPPIDistLD] = TellMePPIInfo(PPINevidenceThreshold,allLDDiseaseGenes,allUniqueGenes);

%-------------------------------------------------------------------------------
% AHBA gene coexpression:
%-------------------------------------------------------------------------------
% For mapped SNPs:
AllenMeanCoexpMapped = TellMeAllenCoexp(allUniqueGenes,allMappedDiseaseGenes);

% For LD SNPs:
AllenMeanCoexpLD = TellMeAllenCoexp(allUniqueGenes,allLDDiseaseGenes);

%===============================================================================
% Assimilate results
%===============================================================================
% Now make a table
gene = allUniqueGenes;
resultsTable = table(gene,numGWASMapped,numLDSNPs,percPPIneigh1Mapped,...
                        percPPIneigh1LD,AllenMeanCoexpMapped,AllenMeanCoexpLD);
                % ,numSNPGenes,numEGenes_LD,numSNPGenes_LD,numLDeGeneseQTL,numLDSNPGeneseQTL); % matchingDrugsString

% Sort by column, then by column, etc. in ordered hierarchy:
resultsTable = sortrows(resultsTable,{'numGWASMapped','numLDSNPs','percPPIneigh1Mapped',... %,'meanPPIDistance'
                'percPPIneigh1LD','AllenMeanCoexpMapped',...
                'AllenMeanCoexpLD'},'descend','MissingPlacement','last');

%-------------------------------------------------------------------------------
% Display just with custom columns
customColumns = {'gene','numGWASMapped','numLDSNPs','percPPIneigh1Mapped',...
                'percPPIneigh1LD','AllenMeanCoexpMapped'}; % ,'matchingDrugsString'
display(resultsTable(1:40,ismember(resultsTable.Properties.VariableNames,customColumns)));

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
