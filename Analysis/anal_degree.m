k = sum(Adj);
kU = unique(k);
fk = arrayfun(@(x)sum(k==x),unique(sum(Adj)));
f = figure('color','w');
bar(unique(k),fk)
