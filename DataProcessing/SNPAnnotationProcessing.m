% Idea is to make a nicer version of the SNPAnnotationTable

LDthreshold = 0.5;

%-------------------------------------------------------------------------------
% IMPORT GWAS-BASED INFORMATION:
%-------------------------------------------------------------------------------
fid = fopen('2_1_SNP_identifier_v2.csv','r');
fprintf(1,'Using new file generated Jul-2018\n')
C = textscan(fid,'%s%u%u%u%u%u%s%u%u%u','Delimiter',',','HeaderLines',1);
fclose(fid);
SNP_id = C{1};
isSZP = logical(C{2});
isADHD = logical(C{3});
isASD = logical(C{4});
isBIP = logical(C{5});
isMDD = logical(C{6});
mappedGene = C{7};
mappedGene(strcmp(mappedGene,'0')) = {''}; % remove '0' -> empty
isGWAS = logical(C{8});
isLD = logical(C{9});
isDiabetes = logical(C{10});
SNPAnnotationTable = table(SNP_id,mappedGene,isGWAS,isLD,isSZP,isADHD,isASD,isBIP,isMDD,isDiabetes);
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% For all unique SNPs, we want a column: LDgenes
%-------------------------------------------------------------------------------
[uniqueSNPs,ia] = unique(SNP_id);
mappedGene_compare = mappedGene(ia);
numUniqueSNPs = length(uniqueSNPs);
fprintf(1,'We have %u unique SNPS\n',length(uniqueSNPs));

numUniqueSNPs = 10;
mappedGenes = cell(numUniqueSNPs,1);
LDgenes = cell(numUniqueSNPs,1);
for i = 1:numUniqueSNPs
    theSNP = uniqueSNPs{i};

    % Get the mapped gene for this SNP:
    mappedGenes{i} = SQL_genesForSNPs(theSNP);

    % Get LD SNPs:
    LD_SNPs = SQL_SNP_LD_SNP(theSNP,LDthreshold);

    % Map each SNP to its set of LD genes
    if isempty(LD_SNPs)
        LDgenes{i} = '';
    else
        LDgenes{i} = SQL_genesForSNPs(LD_SNPs);
    end
end

SNPgeneTable = table(SNP_id,mappedGene);
