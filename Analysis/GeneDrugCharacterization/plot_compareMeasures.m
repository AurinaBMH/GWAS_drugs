function f = plot_compareMeasures(whatDiseases_GWAS, whatMeasures, recalc, whatNull)

if nargin < 3
    recalc = false; 
end
if nargin < 4
    whatNull = 'randomDrugR_all_drugbank'; 
end


if recalc
    [Ptable, measureNames] = compare_optimizedScores(whatDiseases_GWAS, whatMeasures, whatNull);
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
    f = figure('color','w', 'Position', [300, 300, 1200, 650]);
    
    % reorder bars by size
    [pPlot, ix] = sort(Ptable.(whatDiseases_GWAS{i}).Pvals, 'descend');
%     switch whatV
%         case 'horizontal'
%             %ax{i} = subplot(1, numDiseases_GWAS, i);
%         case 'vertical'
%             %ax{i} = subplot(numDiseases_GWAS, 1, i);
%     end
    %hold on;
    pPlot_all(i,:) = pPlot; 
    %b = bar(pPlot);
    hold on; 
    stem(pPlot, 'Marker','none', 'LineStyle',':', 'Color',[.25 .25 .25], 'LineWidth',2)
    
    b = scatter(1:length(ix), pPlot, 500, colors(ix,:),...
              'MarkerFaceColor',[1 1 1],...
              'LineWidth',6); 
    set(gcf, 'renderer', 'painters') 
    xticks(1:numMeasures); 
    xticklabels(Mlabels(ix)); 
    xtickangle(90); 

    %ax{i}.XTick = 1:numMeasures;
    %ax{i}.XTickLabel = Mlabels(ix);
    %ax{i}.XTickLabelRotation = 45;
    
 %   b.CData = colors(ix,:);
%    b.FaceColor = 'flat';
%    b.EdgeColor = [1 1 1]; 
    box off; 
    
    set(gca,'FontSize', 25)
    
    ylabel('-log10(P)')
    title(sprintf('%s', whatDiseases_GWAS{i}))
    
    
    %if strcmp(whatV, 'vertical')
    %    if i==numDiseases_GWAS
    %        xlabel('Measures')
    %    end
    %else
    xlabel('Measures')
    %end
    ylim([0 4]);
    % line for BF corrected value, p=0.01
    if strcmp(whatMeasures, 'allPsych')
        yline(-log10(0.05/27), ':', 'color', [.15 .15 .15], 'LineWidth', 3);
        % 27 measures for psych
    elseif strcmp(whatMeasures, 'allBody')
        yline(-log10(0.05/28), ':', 'color', [.15 .15 .15], 'LineWidth', 3);
        % 28 measures for body
    end
    % line for 0.05
    yline(-log10(0.05/6), ':', 'color', [160,160,160]/255, 'LineWidth', 3);
    
    figureName = sprintf('figures/compareMeasures_%s_%s', whatDiseases_GWAS{i}, whatMeasures);
    print(f,figureName,'-dpng','-r300');
end

%figure; imagesc(pPlot_all); 
%linkaxes([ax{:}],'y');

end