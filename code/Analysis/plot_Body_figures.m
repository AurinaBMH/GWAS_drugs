function plot_Body_figures()

params = SetDefaultParams();
similarityTypes = {'MAGMAdefault', 'PPI_mapped_th600', 'eQTLbrain', 'AlleneQTLbrain'};
whatDiseases_GWAS = {'IBD','RA', 'HF', 'DIABETES'}; %'BIPandSCZ'
whatMeasures = 'allBody';
whatNull = sprintf('randomDrugR_%s_drugbank', params.whatTargets);

%-------------------------------------------------------
% Figure S3: Correspondence across different data processing methods
%-------------------------------------------------------
DOrecalc = false;
plot_compareMeasures(whatDiseases_GWAS, whatMeasures, DOrecalc);

%-------------------------------------------------------
% Additional figure (not shown in the manuscript) - pairwise associations between non-psychiatric disorders and their treatment targets
% (same as Figure S1 for psychiatric disorders)
%-------------------------------------------------------

for s=1:length(similarityTypes)
    
    if contains(similarityTypes{s},'PPI')
        whatProperty = 'percPPIneighbors1';
    else
        if ~contains(similarityTypes{s},'Allen')
            whatProperty = 'P';
        elseif contains(similarityTypes{s},'Allen')
            whatProperty = 'zval';
        end
    end
    
    [rhosALL ,pValsALL] = DistinguishingCharBar(similarityTypes{s},whatProperty, whatNull, 'BF', whatDiseases_GWAS, true, length(similarityTypes), whatMeasures);
    figureName = sprintf('figures_2022/BarChart_body_%s_%s_%s', similarityTypes{s}, whatMeasures, whatNull);
    print(gcf,figureName,'-dpng','-r300');
    
    
end
end

