
disease1 = 'SZP';
disease2 = 'BIP';

isD1 = ismember(eQTLproteinnames,...
            unique(eQTLidentifier.Name(eQTLidentifier.(disease1) & ~eQTLidentifier.Partners)));
isD2 = ismember(eQTLproteinnames,...
            unique(eQTLidentifier.Name(eQTLidentifier.(disease2) & ~eQTLidentifier.Partners)));

isEither = (isD1|isD2);
Adj_filter = Adj(isEither,isEither);
isD1_filter = isD1(isEither);
isD2_filter = isD2(isEither);

[~,reOrder] = sort(sum([isD1_filter,isD2_filter],2),'descend');
Adj_filter = Adj_filter(reOrder,reOrder);
isD1_filter = isD1_filter(reOrder);
isD2_filter = isD2_filter(reOrder);

f = figure('color','w');
imagesc([repmat(isD1_filter,1,10),repmat(isD2_filter,1,10),Adj_filter]);
colormap(gray);
axis('square')


diseases = {'ADHD','ASD','BIP','MDD','SZP'};
numDiseases = length(diseases);
howMuchOverlap = zeros(numDiseases,numDiseases);
for i = 1:numDiseases
    for j = i:numDiseases
        howMuchOverlap(i,j) = quantifyOverlap(diseases{i},diseases{j},...
                            eQTLproteinnames,eQTLidentifier,Adj);
    end
end
f = figure('color','w'); ax = gca;
imagesc(howMuchOverlap)
ax.XTick = 1:numDiseases;
ax.XTickLabel = diseases;
ax.YTick = ax.XTick;
ax.YTickLabel = ax.XTickLabel;
