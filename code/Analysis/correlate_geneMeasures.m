% correlate and group all measures
function [f, Mnames, Mnumbers] = correlate_geneMeasures(disorder, whatMeasures, plotSeparate)
% select_measures - allPsych - relevant to psychiatric
% select_measures - allBody - relevant to non-psychiatric
% select_measures - all - all available
if nargin < 2
    whatMeasures = 'allPsych'; 
end

if nargin<3
    plotSeparate = true; 
end

% load data
fileName = sprintf('DataOutput_2024/resultsTable_%s_BF_2024_all_drugbank.mat', disorder);
load(fileName, 'geneScores');

% get all measures
M = fieldnames(geneScores); 
params = SetDefaultParams();

switch whatMeasures
    case 'allPsych'
        measures = params.whatANNOT_psych';
    case 'allBody'
        measures = params.whatANNOT_body';
    case 'all'
        measures = setdiff(M, {'gene', 'params'}); 
end

measureNames = {'numPPIneighbors1', 'percPPIneighbors1'};

% make a vector of all measures
p=1;
for k=1:length(measures)
    
    if contains(measures{k}, 'Allen')
        MN{p} = {measures{k}, 'zval'};
        p=p+1;
    elseif contains(measures{k}, 'PPI')
        for l=1:length(measureNames)
            MN{p} = {measures{k}, measureNames{l}};
            p=p+1;
        end
    else
        MN{p} = {measures{k}, 'P'};
        p=p+1;
    end
end

% calculate correlations between each pair
r = zeros(length(MN));
p = zeros(length(MN));
for i=1:length(MN)
    for j=i+1:length(MN)
        
        m1 = geneScores.(MN{1,i}{1,1}).(MN{1,i}{1,2});
        m2 = geneScores.(MN{1,j}{1,1}).(MN{1,j}{1,2});
        
        [r(i,j),p(i,j)] = corr(m1, m2, 'type', 'Spearman', 'rows', 'complete');
        
        
    end
end
% reorder
r=r+r';
r(boolean(eye(size(r,1)))) = NaN;

r1 = all(isnan(r),2); 
r2 = all(isnan(r),1); 

r(:,r1) = []; 
r(r2,:) = []; 

MN(r1) = []; 

ord = BF_ClusterReorder(r);
RR = r(ord, ord);
%PP = p(ord, ord);


% combine measure names
Mnames = cell(length(MN),1);

for tt=1:length(MN)
    MN{tt} = strrep(MN{tt},'_',' ');
    UU = join(MN{tt}, '.');
    Mnames{tt} = UU{1};
end
    
Mlabels = give_MeasureLabels(Mnames); 

% reorder based on clustered matrix
Mnumbers = 1:length(Mnames); 
Mnames = Mnames(ord);
Mnumbers = Mnumbers(ord); 
Mlabels = Mlabels(ord); 


% plot
colors = cbrewer('div', 'RdBu', 64);
colors = flipud(colors);
% if plotSeparate
%     f = figure('color','w', 'Position', [100 100 1200 1200]);
% end
% imagesc(RR); axis('equal');
% colormap(colors);
% 
% yticks(1:length(MN));
% yticklabels(Mlabels)
% 
% xticks(1:length(MN));
% xticklabels(Mlabels)
% xtickangle(45);
% 
% caxis([-1 1])
% colorbar
% title(sprintf('%s', disorder))
% set(gca,'FontSize', 22)


[measureColors, ~, cluster] = give_measureColors(Mnames);
[~,V] = unique(cluster); 
measureColors = measureColors(V,:); 
% CLUSTER 1: SNPposition - green; 
% CLUSTER 2: Chromatin - pink
% CLUSTER 3: eQTL - green
% CLUSTER 4: AHBA - grey
% CLUSTER 5: PPIperc - light blue
% CLUSTER 6: PPInum - dark blue
   
% assign clusters
if plotSeparate
    f = figure('color','w', 'Position', [100 100 1000 800]);
end

plotClusteredData(RR,cluster,Mlabels',colors,measureColors);
title(sprintf('%s', disorder))
set(gca,'FontSize', 18)

             

end