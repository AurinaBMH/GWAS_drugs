function f = plot_measureOverview(Pmatrix, T, similarityTypes_label)
% make a bar plot measure to disorder

params = SetDefaultParams();
[~, ord] = intersect(params.whatDiseases_Treatment, T, 'stable'); 
whatDiseases_Treatment_label = params.whatDiseases_Treatment_label(ord); 

f = figure('color','w', 'Position', [300, 300, 1500, 400]);
ax = cell(size(Pmatrix,1),1);
cMapGeneric = BF_getcmap('set4',size(Pmatrix,2),false);

for l=1:size(Pmatrix,1)

    cMapGeneric_n = cMapGeneric;

    Pplot = Pmatrix(l,:);
    [Pbar, ix] = sort(-log10(Pplot), 'descend'); 

    ax{l} = subplot(1,size(Pmatrix,1),l); hold on
   
    b = bar(Pbar);
    for k=1:length(Pplot)
        if Pplot(k)>0.05
            cMapGeneric_n(k,:) = brighten(cMapGeneric(k,:),0.99);
        elseif Pplot(k)<=0.05 && Pplot(k)>=0.05/(length(Pplot))
            cMapGeneric_n(k,:) = brighten(cMapGeneric(k,:),0.85);
        end
    end
    
     ax{l}.XTick = 1:length(whatDiseases_Treatment_label);
     ax{l}.XTickLabel = whatDiseases_Treatment_label(ix);
                
    b.CData = cMapGeneric_n(ix,:);
    b.FaceColor = 'flat';
    
    
    ax{l}.XTickLabelRotation = 45;
    xlabel('Disorder')
    ylabel('-log10(P)')    
    title(sprintf('%s',similarityTypes_label{l}),'interpreter','none')
        
        
    linkaxes([ax{:}],'y');
    set(gca,'FontSize', 14)
end




end