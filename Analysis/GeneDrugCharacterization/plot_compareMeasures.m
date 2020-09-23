function plot_compareMeasures(whatDiseases_GWAS, whatMeasures)

[Ptable, measureNames] = compare_optimizedScores(whatDiseases_GWAS, whatMeasures, false);

numDiseases_GWAS = length(whatDiseases_GWAS);
numMeasures = length(measureNames);

% plot all measures in a single bar chart, order by p-value, colour by
% type: PPI, Allen, eQTL, MAGMAh, Possition, Combined
colors = give_measureColors(measureNames);

% one figure;
f = figure('color','w', 'Position', [300, 300, 800, 1500]);

for i=1:numDiseases_GWAS
    % reorder bars by size
    [pPlot, ix] = sort(Ptable.(whatDiseases_GWAS{i}).Pvals, 'descend');
    ax{i} = subplot(numDiseases_GWAS,1,i); hold on;

    b = bar(pPlot);
    ylabel('-log10(P)')
    xlabel('Measures')
    title(sprintf('%s', whatDiseases_GWAS{i}))
    
    ax{i}.XTick = 1:numMeasures;
    %ax{i}.XTickLabel = measureNames(ix);
    
    b.CData = colors(ix,:);
    b.FaceColor = 'flat';
    
    ax{i}.XTickLabelRotation = 45;
    
end

linkaxes([ax{:}],'y');
set(gca,'FontSize', 14)



end