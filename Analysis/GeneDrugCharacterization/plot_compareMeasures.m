function f = plot_compareMeasures(whatDiseases_GWAS, whatMeasures, recalc, whatV)

if nargin < 3
    recalc = false; 
end

if nargin < 4
    whatV = 'horizontal';
end

if recalc
    [Ptable, measureNames] = compare_optimizedScores(whatDiseases_GWAS, whatMeasures, false);
    fileName = sprintf('DataOutput/Ptable_%s.mat', whatMeasures); 
    save(fileName, 'Ptable', 'measureNames'); 
else
    load(sprintf('Ptable_%s.mat', whatMeasures))
end

numDiseases_GWAS = length(whatDiseases_GWAS);
numMeasures = length(measureNames);
Mlabels = give_MeasureLabels(measureNames); 

% plot all measures in a single bar chart, order by p-value, colour by
% type: PPI, Allen, eQTL, MAGMAh, Possition, Combined
[colors, measureLabels] = give_measureColors(measureNames);
measureNR = 1:length(measureLabels); 

% one figure;
% switch whatV
%     case 'horizontal'
%         f = figure('color','w', 'Position', [300, 300, 2500, 300]);
%     case 'vertical'
%         f = figure('color','w', 'Position', [300, 300, 600, 2500]);
% end

for i=1:numDiseases_GWAS
    
    % plot separate figures; 
    f = figure('color','w', 'Position', [300, 300, 1200, 600]);
    
    % reorder bars by size
    [pPlot, ix] = sort(Ptable.(whatDiseases_GWAS{i}).Pvals, 'descend');
    switch whatV
        case 'horizontal'
            %ax{i} = subplot(1, numDiseases_GWAS, i);
        case 'vertical'
            %ax{i} = subplot(numDiseases_GWAS, 1, i);
    end
    %hold on;

    b = bar(pPlot);
    
    xticks(1:numMeasures); 
    xticklabels(Mlabels(ix)); 
    xtickangle(45); 
    %ax{i}.XTick = 1:numMeasures;
    %ax{i}.XTickLabel = Mlabels(ix);
    %ax{i}.XTickLabelRotation = 45;
    
    b.CData = colors(ix,:);
    b.FaceColor = 'flat';
    b.EdgeColor = [1 1 1]; 
    box off; 
    
    set(gca,'FontSize', 16)
    
    ylabel('-log10(P)')
    title(sprintf('%s', whatDiseases_GWAS{i}))
    
    
    if strcmp(whatV, 'vertical')
        if i==numDiseases_GWAS
            xlabel('Measures')
        end
    else
        xlabel('Measures')
    end
    ylim([0 4]); 
end

%linkaxes([ax{:}],'y');

end