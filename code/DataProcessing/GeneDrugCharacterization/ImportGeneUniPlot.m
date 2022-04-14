function [geneNameHGNC,proteinNameUniprot,allUniqueProteins] = ImportGeneUniProt(allUniqueGenes,PPINGeneNames)
%-------------------------------------------------------------------------------
% Match genes to their proteins to help matching to PPIN names?:
% (Janette provided the file: 'HGNCgene_to_UniprotProtein.txt')
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Read in the data:
%-------------------------------------------------------------------------------
fileNameProtein = 'HGNCgene_to_UniprotProtein-UTF8.csv';
fprintf(1,'Reading gene->protein mappings from %s\n',fileNameProtein);
fid = fopen(fileNameProtein,'r');
S = textscan(fid,'%s%s%s','HeaderLines',1);
fclose(fid);
geneNameHGNC = S{1};
proteinNameUniprot = S{3};
fprintf(1,'Done.\n');

%-------------------------------------------------------------------------------
% Map genes -> proteins
%-------------------------------------------------------------------------------
allUniqueProteins = cell(length(allUniqueGenes),1);
for i = 1:length(allUniqueGenes)
    theMatch = strcmp(geneNameHGNC,allUniqueGenes{i});
    if sum(theMatch)==1
        allUniqueProteins{i} = proteinNameUniprot{theMatch};
    else
        warning('No protein name found for %s',allUniqueGenes{i});
        allUniqueProteins{i} = '';
    end
end

%-------------------------------------------------------------------------------
% Check which gene names (from DrugGene data) exist in the PPI data:
%-------------------------------------------------------------------------------
isPPIProtein = ismember(allUniqueProteins,PPINGeneNames);
isPPIGene = ismember(allUniqueGenes,PPINGeneNames);
isPPIGeneOrProtein = (isPPIProtein | isPPIGene);
fprintf(1,'%u/%u genes (%u/%u from proteins, +%u) from DrugGene data are in the PPI network\n',...
            sum(isPPIGene),length(isPPIGene),sum(isPPIProtein),length(isPPIProtein),...
            sum(isPPIGeneOrProtein)-sum(isPPIGene));

%-------------------------------------------------------------------------------
% List proteins with no matches in PPIN:
%-------------------------------------------------------------------------------
noMatch = find(~isPPIGeneOrProtein);
for i = 1:length(noMatch)
    matchingProtein = allUniqueProteins{noMatch(i)};
    if isempty(matchingProtein)
        matchingProtein = '<no-match>';
    end
    fprintf(1,'No match in PPIN for %s [protein: %s]\n',...
                    allUniqueGenes{noMatch(i)},matchingProtein);
end


end
