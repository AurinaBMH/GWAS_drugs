

normalizeWithinDrugs = true;
whichDiseases = {'ADHD','BIP','SZP','MDD'};
indicatorTable = ImportTreatmentLists(normalizeWithinDrugs)


[~,ia,ib] = intersect({'ADHD','MDD','BIP','SZP'},indicatorTable.Properties.VariableNames,'stable');
indicatorTable = indicatorTable(:,ib);

%-------------------------------------------------------------------------------
% Look at the top X genes in particular:
numTop = 25;
% Order by maximum weight:
maxGene = max(indicatorTable{:,:},[],2);
[~,ix] = sort(maxGene,'descend');
indicatorTable = indicatorTable(ix,:);
data = indicatorTable{1:numTop,:}';

%-------------------------------------------------------------------------------
f = figure('color','w');
imagesc(data)
colormap([1,1,1;flipud(BF_getcmap('spectral',9,false))])
% colormap([1,1,1;jet(12)])
ax = gca;
ax.YTick = 1:size(data,1);
ax.YTickLabel = indicatorTable.Properties.VariableNames;
ax.XTick = 1:numTop;
ax.XTickLabel = indicatorTable.Row(1:numTop);
ax.XTickLabelRotation = 90;
cB = colorbar;
cB.Label.String = 'Gene treatment weight (arbitrary scale)';
title('Drug weightings over genes')
