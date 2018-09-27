% Parameters:
numTop = 30; % Look at the top X genes in particular
normalizeWithinDrugs = true;
whichDiseases = {'ADHD','MDD','BIP','SZP'};
% 'diabetes','pulmonary','cardiology','gastro'};

%-------------------------------------------------------------------------------
indicatorTable = ImportTreatmentLists(normalizeWithinDrugs)
[~,ia,ib] = intersect(whichDiseases,indicatorTable.Properties.VariableNames,'stable');
indicatorTable = indicatorTable(:,ib);

%-------------------------------------------------------------------------------
% Normalize:
dataNorm = indicatorTable{:,:}';
for i = 1:4
    dataNorm(i,:) = dataNorm(i,:)/norm(dataNorm(i,:),2);
end
% Order by maximum weight:
maxGene = max(dataNorm,[],1);
[~,ix] = sort(maxGene,'descend');
data = dataNorm(:,ix(1:numTop));

%-------------------------------------------------------------------------------
f = figure('color','w');
imagesc(data)
colormap([1,1,1;BF_getcmap('blues',9,false)])
% colormap([1,1,1;jet(12)])
ax = gca;
ax.YTick = 1:size(data,1);
ax.YTickLabel = indicatorTable.Properties.VariableNames;
ax.XTick = 1:numTop;
ax.XTickLabel = indicatorTable.Row(1:numTop);
ax.XTickLabelRotation = 90;
cB = colorbar;
cB.Label.String = 'Normalized treatment weight';
caxis([0,0.4])
title('Drug weightings over genes')
