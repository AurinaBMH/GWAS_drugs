%---------------------------------------------------------------------------
% This script imports HMAGMA output files and assigns each gene an entrezID
%---------------------------------------------------------------------------
clear all;
close all;

params = SetDefaultParams();
Disorders = params.whatGWAS; 
whatANNOT = params.whatANNOT; 

% MAGMA: 'MAGMAdefault'
% MAGMA-H: 'Adult_brain', 'Fetal_brain', 'Neuro', 'Astro',
% eMAGMA psychENCODE: 'eQTLbrain'; 
% eMAGMA GTEx: 'Small_Intestine_Terminal_Ileum' 'Pancreas' 'Whole_Blood' 'Liver' 'Heart_Left_Ventricle' 'Colon_Transverse' 'Colon_Sigmoid' 'Adipose_Subcutaneous' 'Adipose_Visceral_Omentum'

% load gene ID matching file and select genes that have entrezIDs
entrezID = importGENEIDfile('data/GWASlists/BIOMART_geneIDs.txt');
entrezID = entrezID(~isnan(entrezID.entrezgene_id),:);

entrezIDonly = entrezID(:, {'ensembl_gene_id', 'entrezgene_id', 'external_gene_name'});

% there are cases where different gene IDs have the same entrez ID and gene name, selet only usique instances
% of entrezID and gene name combinations and - these are the IDs that will
% be used later.
entrezIDname = entrezID(:, {'entrezgene_id', 'external_gene_name'});
[~, INDunique] = unique(entrezIDname, 'rows');
entrezIDonly = entrezIDonly(INDunique,:);

% make a structure to save the formated outputs;
DISORDERlist=struct;

for D=1:length(Disorders)
    for A=1:length(whatANNOT)

        fileName = sprintf('data/GWASlists/GWASgenes/pgc%s_%s_genes.txt', Disorders{D}, whatANNOT{A});

        %GTEx - based files have genes labeled with entrezIDs, other with stable IDs
        isGTEx = strcmp(whatANNOT{A}, 'eQTLWhole_Blood') || strcmp(whatANNOT{A}, 'eQTLLiver') || strcmp(whatANNOT{A}, 'eQTLHeart_Left_Ventricle') || ...
            strcmp(whatANNOT{A}, 'eQTLPancreas') || strcmp(whatANNOT{A}, 'eQTLSmall_Intestine_Terminal_Ileum') || strcmp(whatANNOT{A}, 'eQTLColon_Transverse') || ...
            strcmp(whatANNOT{A}, 'eQTLColon_Sigmoid') || strcmp(whatANNOT{A}, 'eQTLAdipose_Subcutaneous') || strcmp(whatANNOT{A}, 'eQTLAdipose_Visceral_Omentum');

        if isGTEx
            mapLIST = importeMAGMAGTExfile(fileName);
            [~, INDmap, indENT] = intersect(mapLIST.GENE, entrezIDonly.entrezgene_id, 'stable');
        else
            mapLIST = importHMAGMAoutfile(fileName);
            [~, INDmap, indENT] = intersect(mapLIST.GENE, entrezIDonly.ensembl_gene_id, 'stable');
        end

        fprintf('Processing disorder: %s, annotation %s\n', Disorders{D}, whatANNOT{A});

        % select only overlapping genes from the table
        % in some cases one gene has multiple entrez IDs, probably related
        % to recent updates, hard to tell which one to choose, so we'll
        % keep gene names as well;
        mapLIST = mapLIST(INDmap,:);

        if isGTEx
            % for non-GTEx, give entrezID the contents of gene and add stable IDs
            ENTREZID = mapLIST.GENE;
            mapLIST.GENE = entrezIDonly.ensembl_gene_id(indENT);
        else
            ENTREZID = entrezIDonly.entrezgene_id(indENT);
        end

        GENENAME = entrezIDonly.external_gene_name(indENT);
        mapLIST = addvars(mapLIST,ENTREZID,GENENAME,'After','GENE');


        % divided
        pThr = 0.05/size(mapLIST,1);
        fprintf('%d genes\n', sum(mapLIST.P<pThr));

        DISORDERlist.(whatANNOT{A}).(Disorders{D}) = mapLIST;
        % there are also some cases where the same ENTREZis and same gene
        % name corresponds to several GENEIDs, if that's the case, keep one

        writetable(mapLIST, sprintf('data/GWASlists/GWASgenes/pgc%s_%s_genes_entrezID.txt', ...
            Disorders{D}, whatANNOT{A}), 'Delimiter','\t');
    end
end

% save sigle file for future analyses
save('DataOutput_2022/GWAS_disordersMAGMA.mat', 'DISORDERlist')
