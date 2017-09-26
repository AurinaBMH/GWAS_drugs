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
fprintf(1,'%u unique genes with drug targets (in Drugbank)\n',numUniqueGenes);
numUniqueDrugs = length(unique(geneDrug.drugName));
fprintf(1,'%u unique drugs with gene targets (in Drugbank)\n',numUniqueDrugs);

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

%-------------------------------------------------------------------------------
% Generate a SNP->Gene map
%-------------------------------------------------------------------------------
SNPGeneMap = SNPAnnotationTable(cellfun(@(x)~isempty(x),SNPAnnotationTable.mappedGeneSNP),1:2);
% Make unique:
[~,ix] = unique(SNPGeneMap.SNP_id);
SNPGeneMap = SNPGeneMap(ix,:);

%-------------------------------------------------------------------------------
% Infer the LD gene (i.e., the gene causing the annotation) for LD annotations
% in SNPAnnotationTable
%-------------------------------------------------------------------------------
% Load the LD relationship data:
fid = fopen('2_2_LD_r2.csv','r');
C = textscan(fid,'%s%f%s','Delimiter',',','HeaderLines',1);
fclose(fid);
SNP_id_1 = C{1};
SNP_id_2 = C{3};
LD = C{2};
LDRelateTable = table(SNP_id_1,SNP_id_2,LD);

% For LD SNPs, we can annotate the SNP that caused the annotation
LDGene = cell(numAnnotations,1);
for i = 1:numAnnotations
    if SNPAnnotationTable.isLD(i)
        % We want to determine which gene caused the annotation:
        theLDSNP1 = LDRelateTable.SNP_id_2(strcmp(LDRelateTable.SNP_id_1,SNPAnnotationTable.SNP_id{i}));
        theLDSNP2 = LDRelateTable.SNP_id_1(strcmp(LDRelateTable.SNP_id_2,SNPAnnotationTable.SNP_id{i}));
        theLDSNPs = vertcat(theLDSNP1{:},theLDSNP2{:});
        if ischar(theLDSNPs), theLDSNPs = {theLDSNPs}; end

        % Match to genes:
        theLDGenes = cellfun(@(x)SNPGeneMap.mappedGeneSNP(strcmp(SNPGeneMap.SNP_id,x)),theLDSNPs,'UniformOutput',false);

        if isempty(theLDGenes)
            warning('No LD gene for %s',SNPAnnotationTable.SNP_id{i})
            LDGene{i} = '';
        elseif
            
        end

    else
        LDGene{i} = '';
    end
end

%-------------------------------------------------------------------------------
% Now loop through genes to characterize...
numGWASMapped = zeros(numUniqueGenes,1);

for i = 1:numUniqueGenes
    fprintf(1,'[%u/%u]: %s\n',i,numUniqueGenes,allUniqueGenes{i});
    % --numGWASMapped: the number of GWAS SNPs mapped directly to a gene
    theSNPs = SNPAnnotationTable.SNP_id(strcmp(SNPAnnotationTable.mappedGeneSNP,allUniqueGenes{i}));
    if isempty(theSNPs)
        numGWASMapped(i) = 0;
    else
        numGWASMapped(i) = length(unique(theSNPs));
    end

    %-------------------------------------------------------------------------------

end
