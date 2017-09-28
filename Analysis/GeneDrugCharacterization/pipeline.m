% Pipeline for producing table characterizing individual genes

%===============================================================================
% LOAD AND PROCESS DATA
%===============================================================================

%-------------------------------------------------------------------------------
% Get genes
%-------------------------------------------------------------------------------
% List of all genes with drug targets in Drugbank are in **gene_ATC_matrix.csv**

fid = fopen('8_1_drug_gene.csv','r');
C = textscan(fid,'%s%s','Delimiter',',');
fclose(fid);
geneName = C{1}; drugName = C{2};
geneDrug = table(geneName,drugName);
allUniqueGenes = unique(geneDrug.geneName);
numUniqueGenes = length(allUniqueGenes);
fprintf(1,'%u genes have drug targets in Drugbank\n',numUniqueGenes);
numUniqueDrugs = length(unique(geneDrug.drugName));
fprintf(1,'%u drugs have gene targets in Drugbank\n',numUniqueDrugs);

%-------------------------------------------------------------------------------
% Read in SNP annotations
%-------------------------------------------------------------------------------
fid = fopen('2_1_SNP_identifier.csv','r');
C = textscan(fid,'%s%u%u%u%u%u%s%u%u','Delimiter',',','HeaderLines',1);
fclose(fid);
SNP_id = C{1};
isSZP = logical(C{2});
isADHD = logical(C{3});
isASD = logical(C{4});
isBIP = logical(C{5});
isMDD = logical(C{6});
mappedGeneSNP = C{7};
mappedGeneSNP(strcmp(mappedGeneSNP,'0')) = {''}; % remove '0' -> empty
isGWAS = logical(C{8});
isLD = logical(C{9});
SNPAnnotationTable = table(SNP_id,mappedGeneSNP,isGWAS,isLD,isSZP,isADHD,isASD,isBIP,isMDD);
numAnnotations = height(SNPAnnotationTable);
fprintf(1,'%u annotations for %u GWAS mapped and %u LD\n',numAnnotations,sum(isGWAS),sum(isLD));
fprintf(1,'%u/%u genes with drug targets have annotations\n',...
                sum(ismember(allUniqueGenes,SNPAnnotationTable.mappedGeneSNP)),numUniqueGenes);

hasGeneName = cellfun(@(x)~isempty(x),mappedGeneSNP);
fprintf(1,'%u/%u annotations are mapped to gene names\n',sum(hasGeneName),numAnnotations);
fprintf(1,'%u/%u GWAS-mapped annotations have gene names\n',sum(isGWAS & hasGeneName),sum(isGWAS));
fprintf(1,'%u/%u LD annotations have gene names\n',sum(isLD & hasGeneName),sum(isLD));

%-------------------------------------------------------------------------------
% Generate a SNP->Gene map
%-------------------------------------------------------------------------------
SNPGeneMap = SNPAnnotationTable(cellfun(@(x)~isempty(x),SNPAnnotationTable.mappedGeneSNP),1:2);
% Make unique:
[~,ix] = unique(SNPGeneMap.SNP_id);
SNPGeneMap = SNPGeneMap(ix,:);
fprintf(1,'%u/%u SNPs have a matching gene name (%u genes)\n',height(SNPGeneMap),...
        length(unique(SNPAnnotationTable.SNP_id)),length(unique(SNPGeneMap.mappedGeneSNP)));

%-------------------------------------------------------------------------------
% Load the LD relationship data:
%-------------------------------------------------------------------------------
fid = fopen('2_2_LD_r2.csv','r');
C = textscan(fid,'%s%f%s','Delimiter',',','HeaderLines',1);
fclose(fid);
SNP_id_1 = C{1};
SNP_id_2 = C{3};
LD = C{2};
LDRelateTable = table(SNP_id_1,SNP_id_2,LD);
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Get cis-eQTL data
%-------------------------------------------------------------------------------
fid = fopen('4_GWAS_eQTL_combn.csv','r');
C = textscan(fid,'%s%s%f%s%f','Delimiter',',','HeaderLines',1);
fclose(fid);
eQTL = C{1};
SNP = C{2};
eGene_qVal = C{3};
eGene_qVal(eGene_qVal==-999) = NaN;
tissueSample = C{4};
SNPgene_pVal = C{5};
SNPgene_pVal(SNPgene_pVal==-999) = NaN;
eQTLTable = table(eQTL,SNP,eGene_qVal,tissueSample,SNPgene_pVal);
isMissing = strcmp(eQTL,'0');
fprintf(1,'%u eQTLs assigned ''0'' -- no hgnc gene symbol -- removed\n',sum(isMissing));
eQTLTable = eQTLTable(~isMissing,:);
fprintf(1,'%u eQTL annotations loaded\n',height(eQTLTable));
fprintf(1,'%u eGene annotations\n',sum(~isnan(eQTLTable.eGene_qVal)));
fprintf(1,'%u SNPgene annotations\n',sum(~isnan(eQTLTable.SNPgene_pVal)));

%-------------------------------------------------------------------------------
% Infer the LD gene (i.e., the gene causing the annotation) for LD annotations
% in SNPAnnotationTable
%-------------------------------------------------------------------------------

% For LD SNPs, we can annotate the SNP that caused the annotation
% LDGene = cell(numAnnotations,1);
% for i = 1:numAnnotations
%     if SNPAnnotationTable.isLD(i)
%         % We want to determine which gene caused the annotation:
%         theLDSNP1 = LDRelateTable.SNP_id_2(strcmp(LDRelateTable.SNP_id_1,SNPAnnotationTable.SNP_id{i}));
%         theLDSNP2 = LDRelateTable.SNP_id_1(strcmp(LDRelateTable.SNP_id_2,SNPAnnotationTable.SNP_id{i}));
%         theLDSNPs = unique(vertcat(theLDSNP1,theLDSNP2));
%         if ischar(theLDSNPs), theLDSNPs = {theLDSNPs}; end
%
%         % Match to genes:
%         theLDGenes = cellfun(@(x)SNPGeneMap.mappedGeneSNP(strcmp(SNPGeneMap.SNP_id,x)),theLDSNPs,'UniformOutput',false);
%
%         if isempty(theLDGenes)
%             % No matches? Must be an LD relationship to a SNP without a gene
%             warning('No LD gene for %s',SNPAnnotationTable.SNP_id{i});
%             LDGene{i} = '';
%         elseif length(theLDGenes)==1
%             % One match -- easy:
%             LDGene{i} = theLDGenes{1};
%         else
%             fprintf(1,'Multiple matches :-O\n');
%         end
%
%     else
%         LDGene{i} = '';
%     end
% end

%-------------------------------------------------------------------------------
% Now loop through genes to characterize...
numGWASMapped = zeros(numUniqueGenes,1);
numLDSNPs = zeros(numUniqueGenes,1);
numEGenes = zeros(numUniqueGenes,1);
numSNPGenes = zeros(numUniqueGenes,1);

for i = 1:numUniqueGenes
    gene_i = allUniqueGenes{i};
    fprintf(1,'[%u/%u]: %s\n',i,numUniqueGenes,gene_i);
    %-------------------------------------------------------------------------------
    % --numGWASMapped: the number of GWAS SNPs mapped directly to a gene
    %-------------------------------------------------------------------------------
    theSNPs = SNPAnnotationTable.SNP_id(strcmp(SNPAnnotationTable.mappedGeneSNP,gene_i) & SNPAnnotationTable.isGWAS);
    if isempty(theSNPs)
        numGWASMapped(i) = 0;
    else
        numGWASMapped(i) = length(unique(theSNPs));
    end

    %-------------------------------------------------------------------------------
    % --numLD: the number of GWAS SNPs LD to a gene
    %-------------------------------------------------------------------------------
    % If I want to count the number of SNPs that are LD to a given gene, I need to:
    % (i) Check every SNP in that gene (?? Using SNP_identifier.csv??)
    % (ii) Get a list of all (unique) SNPs are LD to these SNPs (Using LD_r2.csv)
    % (iii) Check all LD SNPs for direct disease annotations (in SNP_identifier.csv).
    isGeneAndLD = strcmp(SNPAnnotationTable.mappedGeneSNP,gene_i) & SNPAnnotationTable.isLD;
    % NB: No way to double check that all LD annotations are from unique SNPs...
    % (is an assumption that Janette has formatted the file in this way)
    numLDSNPs(i) = sum(isGeneAndLD);

    %-------------------------------------------------------------------------------
    % eQTLs-GWAS
    %-------------------------------------------------------------------------------
    % eGenes
    numEGenes(i) = sum(~isnan(eQTLTable.eGene_qVal) & strcmp(eQTLTable.eQTL,gene_i));
    % SNPGenes
    numSNPGenes(i) = sum(~isnan(eQTLTable.SNPgene_pVal) & strcmp(eQTLTable.eQTL,gene_i));

    %-------------------------------------------------------------------------------
    % eQTLs-LD
    %-------------------------------------------------------------------------------

end

% Now make a table
gene = allUniqueGenes;
results = table(gene,numGWASMapped,numLDSNPs,numEGenes,numSNPGenes);
results = sortrows(results,{'numGWASMapped','numLDSNPs','numEGenes','numSNPGenes'},'descend');
fprintf(1,'%u genes have mapped and LD\n',sum(results.numLDSNPs>0 & results.numGWASMapped>0));
fprintf(1,'%u genes have mapped but no LD\n',sum(results.numLDSNPs==0 & results.numGWASMapped>0));
fprintf(1,'%u genes have LD but no mapped\n',sum(results.numLDSNPs>0 & results.numGWASMapped==0));
fprintf(1,'%u genes have eGene annotations\n',sum(results.numEGenes>0));
fprintf(1,'%u genes have SNPgene annotations\n',sum(results.numSNPGenes>0));
