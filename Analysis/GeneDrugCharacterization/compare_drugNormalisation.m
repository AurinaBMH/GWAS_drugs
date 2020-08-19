% look for correlations between drug-normalised and non-normalised scores;
function compare_drugNormalisation()

[indicatorTableF] = ImportTreatmentLists(false);
[indicatorTableT] = ImportTreatmentLists(true);
geneNamesF = indicatorTableF.Row;
geneNamesT = indicatorTableT.Row;
[~, indF, indT] = intersect(geneNamesF, geneNamesT);

indicatorTableF = indicatorTableF(indF,:);
indicatorTableT = indicatorTableT(indT,:);

disorders = indicatorTableF.Properties.VariableNames;

figure('color', 'w');
for i=1:size(indicatorTableF,2)
    
    [r,p] = corr(indicatorTableF.(disorders{i}), indicatorTableT.(disorders{i}));
    subplot(2,4,i)
    scatter(indicatorTableF.(disorders{i}), indicatorTableT.(disorders{i}), ...
        40, 'MarkerEdgeColor',[0 .5 .5],...
        'MarkerFaceColor',[0 .7 .7],...
        'LineWidth',1.5)
    title(sprintf('%s, r=%.2f, p=%.2d', disorders{i}, r, p))
    xlabel('non-normalised scores')
    ylabel('drug-normalised scores')
    
end

end

