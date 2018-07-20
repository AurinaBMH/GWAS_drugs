function resultsTable = pipeline(whatDisease,PPINevidenceThreshold,LDThreshold)
% Pipeline for producing table characterizing individual genes
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Parse inputs:
%-------------------------------------------------------------------------------
if nargin < 1
    % GWAS hits from which disease?: 'all','SZP','ASD','ADHD','BIP','MDD'
    whatDisease = 'all';
end
if nargin < 2
    PPINevidenceThreshold = 0; % evidence threshold for including PPI interactions
end
if nargin < 3
    LDThreshold = 0.5;
end

%===============================================================================
% LOAD DATA
%===============================================================================

%-------------------------------------------------------------------------------
% Restrict our characterization to genes associated with (any) drug-treatment:
%-------------------------------------------------------------------------------
% List of all genes with drug targets in Drugbank are in **gene_ATC_matrix.csv**
% [geneDrugTable,allUniqueGenes] = DrugGeneImport();
% Import drug classification
% drugClassTable = DrugClassImport();
indicatorTable = ImportTreatmentLists(false);
allUniqueGenes = indicatorTable.Properties.RowNames;
numUniqueGenes = length(allUniqueGenes);

%-------------------------------------------------------------------------------
% Import SNP, gene, disease, GWAS, LD annotations for a given GWAS study:
%-------------------------------------------------------------------------------
% [SNPAnnotationTable,SNPGeneMap,allDiseaseSNPs] = SNPAnnotationImport(whatDisease);
fileName = sprintf('SNPAnnotationTable_%u.mat',LDthreshold*100);
load(fileName,'SNPAnnotationTable');

% Combine all mapped genes:
allMappedDiseaseGenes = unique(vertcat(SNPAnnotationTable.mappedGenes{:}));

% Get LD-genes:
% isLDAndNotEmpty = SNPAnnotationTable.isLD & cellfun(@(x)~isempty(x),SNPAnnotationTable.mappedGene);
% allDiseaseGenesLD = unique(SNPAnnotationTable.mappedGene(isLDAndNotEmpty));

% fprintf(1,'%u/%u genes with drug targets have annotations\n',...
%         sum(ismember(allUniqueGenes,SNPAnnotationTable.mappedGene)),numUniqueGenes);

%===============================================================================
%======================CHARACTERIZATION MODULES=================================
%===============================================================================

% Now get different characterizations of each gene of interest to closeness
% to different data types:
[numGWASMapped,numLDSNPs] = TellMeDNADistance(allUniqueGenes,SNPAnnotationTable);

%-------------------------------------------------------------------------------
% eQTL Characterization:
% eQTL_info = TellMeQTL(allDiseaseGenesMapped,allUniqueGenes);
% numEGenes = eQTL_info.numEGenes;
% numSNPGenes = eQTL_info.numSNPGenes;
% numEGenes_LD = eQTL_info.numEGenes_LD;
% numSNPGenes_LD = eQTL_info.numSNPGenes_LD;
% numLDeGeneseQTL = eQTL_info.numLDeGeneseQTL;
% numLDSNPGeneseQTL = eQTL_info.numLDSNPGeneseQTL;

%-------------------------------------------------------------------------------
% PPIN characterization:
%-------------------------------------------------------------------------------

% Just mapped disease genes (genes with a GWAS SNP in them):
[numPPIneigh1Mapped,percPPIneigh1Mapped,meanPPIDistMapped] = TellMePPIInfo(PPINevidenceThreshold,allDiseaseGenesMapped,allUniqueGenes);

% Just genes LD to GWAS SNPs:
[numPPIneigh1LD,percPPIneigh1LD,meanPPIDistLD] = TellMePPIInfo(PPINevidenceThreshold,allDiseaseGenesLD,allUniqueGenes);

%-------------------------------------------------------------------------------
% AHBA gene coexpression:
%-------------------------------------------------------------------------------
% For mapped SNPs:
AllenMeanCoexpMapped = TellMeAllenCoexp(allUniqueGenes,allDiseaseGenesMapped);

% For LD SNPs:
AllenMeanCoexpLD = TellMeAllenCoexp(allUniqueGenes,allDiseaseGenesLD);

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
                % 'AllenMeanCoexp',...
                % 'numEGenes','numSNPGenes','numEGenes_LD',...
                % 'numSNPGenes_LD','numLDeGeneseQTL','numLDSNPGeneseQTL'}


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
