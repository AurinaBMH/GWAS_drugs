% convert gene IDs for enrichment and save to file
function save_enrichment_scores()

entrezIDs = readtable('enrichment_2024/SynGO_id_convert_2024-02-04_1415/idmap.xlsx');
similarityTypes = {'PPI_mapped_th600'}; %{'MAGMAdefault', 'PPI_mapped_th600', 'eQTLbrain', 'AlleneQTLbrain'};
disorders = {'ADHD', 'BIP', 'MDD', 'SCZ', 'DIABETES'};

for i=1:length(similarityTypes)
    for d=1:length(disorders)
        
        load(sprintf('enrichment_2024/enrichment_drug_%s', similarityTypes{i}));
        load(sprintf('enrichment_2024/enrichment_GWAS_%s', similarityTypes{i}));
        
        [~, Ie, Idrug] = intersect(entrezIDs.query, enrichment_score_drug.geneName);
        
        % for drugs
        enrichment_score_drug.geneName(Idrug) = entrezIDs.entrezgene(Ie);
        
        drug_names = enrichment_score_drug.Properties.VariableNames;
        
        TF1 = contains(drug_names,disorders{d});
        TF2 = contains(drug_names,'geneName');
        TF = TF1|TF2;
        
        enrichment_score_drug_dis = enrichment_score_drug(:,TF);
        writetable(enrichment_score_drug_dis, sprintf('enrichment_2024/enrichment_%s_drug_%s', disorders{d}, similarityTypes{i}), 'Delimiter', '\t')
        
        % for GWAS
        gwas_names = fieldnames(enrichment_score_GWAS);
        GT = contains(gwas_names, disorders{d});
        
        enrichment_score_gwas_dis = enrichment_score_GWAS.(gwas_names{GT});
        [~, Ie, Igwas] = intersect(entrezIDs.query, enrichment_score_gwas_dis.geneNames);
        
        enrichment_score_gwas_dis.geneNames(Igwas) = entrezIDs.entrezgene(Ie);
        writetable(enrichment_score_gwas_dis, sprintf('enrichment_2024/enrichment_%s_gwas_%s', disorders{d}, similarityTypes{i}), 'Delimiter', '\t')
        
    end
end
end

