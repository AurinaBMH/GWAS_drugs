% correlate and group all measures
function correlate_geneMeasures(disorder)
if nargin<1
    disorder = 'DIABETES';
end

fileName = sprintf('resultsTable_%s_BF_2020.mat', disorder);
load(fileName, 'geneScores');

% get all measures
measures = setdiff(fieldnames(geneScores), {'params', 'gene'});
measureNames = {'numPPIneighbors1', 'percPPIneighbors1', 'numPPIneighbors2', 'percPPIneighbors2'};

% make a vector of all measures
p=1;
for k=1:length(measures)
    
    if contains(measures{k}, 'Allen')
        MN{p} = {measures{k}};
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

r = zeros(length(MN));
p = zeros(length(MN));
for i=1:length(MN)
    for j=i+1:length(MN)
        
        
        if contains(MN{1,i}{1,1}, 'Allen')
            m1 = geneScores.(MN{1,i}{1,1});
        else
            m1 = geneScores.(MN{1,i}{1,1}).(MN{1,i}{1,2});
        end
        
        if contains(MN{1,j}{1,1}, 'Allen')
            m2 = geneScores.(MN{1,j}{1,1});
        else
            m2 = geneScores.(MN{1,j}{1,1}).(MN{1,j}{1,2});
        end
        
        [r(i,j),p(i,j)] = corr(m1, m2, 'type', 'Spearman', 'rows', 'complete');
        
        
    end
end

r=r+r';
r(boolean(eye(size(r,1)))) = NaN;
ord = BF_ClusterReorder(r);
RR = r(ord, ord);
PP = p(ord, ord);


% combine measure names
Mnames = cell(length(MN),1);
for tt=1:length(MN)
    MN{tt} = strrep(MN{tt},'_',' ');
    UU = join(MN{tt}, '.');
    Mnames{tt} = UU{1};
end
% reorder based on clustered matrix
Mnames = Mnames(ord);

colors = cbrewer('div', 'RdBu', 64);
colors = flipud(colors);
figure('color','w', 'Position', [100 100 1200 1200]);
imagesc(RR); axis('square')
colormap(colors);

yticks(1:length(MN));
yticklabels(Mnames)

xticks(1:length(MN));
xticklabels(Mnames)
xtickangle(45);

caxis([-1 1])
colorbar
title(sprintf('%s', disorder))


end