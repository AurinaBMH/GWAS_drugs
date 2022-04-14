%----------------------------------------------------------------------
% This script loads raw protein link data and replaces protein IDs with
% corresponding gene names based on biomart annotation.
%----------------------------------------------------------------------
function make_PPI_linkfile()

% load PPI network links
PROTEINlinks = importPROTEINLINKSfile('9606.protein.links.v11.0.txt');
P = PROTEINlinks(:,[1,2]);
nPairsR = unique(P);
if size(nPairsR,1)==size(PROTEINlinks,1)
    sprintf('all protein pairs are unique\n')
else
    sprintf('there are duplicate protin pairs')
end
%  - at the level of protein all pairs are unique
% remove 9606 from protein IDs and replace protein IDs
splitP1 = split(PROTEINlinks.protein1,'.');
splitP2 = split(PROTEINlinks.protein2,'.');

PROTEINlinks.protein1 = splitP1(:,2);
PROTEINlinks.protein2 = splitP2(:,2);

% load biomart file containing gene names associated with proteins
BioMartRESULTS = importPROTEINfile('data/GWASlists/BIOMART_geneIDs.txt');

% find gene names for proteins and create a table for protein-protein
% interractions with gene names

prot = unique(vertcat(PROTEINlinks.protein1,PROTEINlinks.protein2));
PPI = PROTEINlinks;
PPI.protein1(:) = '';
PPI.protein2(:) = '';

for p=1:length(prot)
    fprintf('Replacing protein nr %d: %s\n', p, prot(p))
    % find where each protein was mentioned
    
    IP1 = strcmp(PROTEINlinks.protein1 ,prot(p));
    IP2 = strcmp(PROTEINlinks.protein2 ,prot(p));
    % replace protein name with corresponding gene name
    IG = strcmp(BioMartRESULTS.ensembl_peptide_id ,prot(p));
    
    if isempty(find(IG, 1))
        PPI.protein1(IP1) = 'uncharacterised';
        PPI.protein2(IP2) = 'uncharacterised';
    else
        PPI.protein1(IP1) = BioMartRESULTS.external_gene_name{IG};
        PPI.protein2(IP2) = BioMartRESULTS.external_gene_name{IG};
    end
end

writetable(PPI,'data/PPIdata/PPIlinks_v11.0.txt', 'WriteVariableNames', true)
end


