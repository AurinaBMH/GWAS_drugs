function pPlot_all = plot_compareMeasures(whatDiseases_GWAS, whatMeasures, recalc, whatYear, whatNull)

if nargin < 3
    recalc = false; 
end

if nargin < 4
    whatYear = '2024'; 
end

if nargin < 5
    whatNull = 'randomDrugR_all_drugbank'; 
end



if recalc
    [Ptable, measureNames] = compare_optimizedScores(whatDiseases_GWAS, whatMeasures, whatNull);
    fileName = sprintf('DataOutput_2024/Ptable_%s_%s.mat', whatMeasures, whatYear); 
    save(fileName, 'Ptable', 'measureNames'); 
else
    load(sprintf('Ptable_%s_%s.mat', whatMeasures, whatYear))
end

numDiseases_GWAS = length(whatDiseases_GWAS);
numMeasures = length(measureNames);
Mlabels = give_MeasureLabels(measureNames); 

% plot all measures in a single bar chart, order by p-value, colour by
% type: PPI, Allen, eQTL, MAGMAh, Possition, Combined
[colors, measureLabels] = give_measureColors(measureNames);
measureNR = 1:length(measureLabels); 

for i=1:numDiseases_GWAS
    
    % plot separate figures; 
    f = figure('color','w', 'Position', [300, 300, 1200, 650]);
    
    % reorder bars by size
    [pPlot, ix] = sort(Ptable.(whatDiseases_GWAS{i}).Pvals, 'descend');
    pPlot_all(i,:) = pPlot; 
    hold on; 
    stem(pPlot, 'Marker','none', 'LineStyle',':', 'Color',[.25 .25 .25], 'LineWidth',2)
    
    b = scatter(1:length(ix), pPlot, 500, colors(ix,:),...
              'MarkerFaceColor',[1 1 1],...
              'LineWidth',6); 
    set(gcf, 'renderer', 'painters') 
    xticks(1:numMeasures); 
    xticklabels(Mlabels(ix)); 
    xtickangle(90); 

    box off; 
    
    set(gca,'FontSize', 25)
    
    ylabel('-log10(P)')
    title(sprintf('%s', whatDiseases_GWAS{i}))
    
    xlabel('Measures')
    %end
    ylim([0 4]);
    % line for BF corrected value
    if strcmp(whatMeasures, 'allPsych')
        yline(-log10(0.05/27), ':', 'color', [.15 .15 .15], 'LineWidth', 3);
        % 27 measures for psych
        yline(-log10(0.05/(27*5)), ':', 'color', [.05 .05 .05], 'LineWidth', 3);
        % 27 measures for psych x 5 disorders
        yline(-log10(0.05/6), ':', 'color', [160,160,160]/255, 'LineWidth', 3);
        % line for 6 measure types
    elseif strcmp(whatMeasures, 'allBody')
        yline(-log10(0.05/28), ':', 'color', [.15 .15 .15], 'LineWidth', 3);
        % 28 measures for body
        yline(-log10(0.05/(28*4)), ':', 'color', [.05 .05 .05], 'LineWidth', 3);
        % 28 measures for body x 4 disorders
        yline(-log10(0.05/5), ':', 'color', [160,160,160]/255, 'LineWidth', 3);
        % line for 5 measure types
    end
    
    figureName = sprintf('figures_%s/compareMeasures_%s_%s', whatYear, whatDiseases_GWAS{i}, whatMeasures);
    print(f,figureName,'-dpng','-r300');
end


end