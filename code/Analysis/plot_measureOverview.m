function f = plot_measureOverview(Pmatrix, T, similarityTypes_label)
% make a bar plot measure to disorder

params = SetDefaultParams();
[~, ord] = intersect(params.whatDiseases_Treatment, T, 'stable'); 
whatDiseases_Treatment_label = params.whatDiseases_Treatment_label(ord); 

for kk=1:length(whatDiseases_Treatment_label)
    whatDiseases_Treatment_label{kk} = strcat(whatDiseases_Treatment_label{kk}, "    "); 
end

f = figure('color','w', 'Position', [300, 300, 1500, 500]);
ax = cell(size(Pmatrix,1),1);
cMapGeneric = BF_getcmap('set4',size(Pmatrix,2),false);

for l=1:size(Pmatrix,1)

    cMapGeneric_n = cMapGeneric;

    Pplot = Pmatrix(l,:);
    Pplot(Pplot==0) = 1/params.numNull; 
    %[Pbar, ix] = sort(-log10(Pplot), 'descend'); 
    % don't order
    Pbar = -log10(Pplot); 
    ix = 1:length(Pplot); 

    ax{l} = subplot(1,size(Pmatrix,1),l); hold on
    % line for BF corrected value, p=0.01
    yline(-log10(0.01), ':', 'color', [.15 .15 .15], 'LineWidth', 3);
    % line for 0.05
    yline(-log10(0.05), ':', 'color', [160,160,160]/255, 'LineWidth', 3);
    %cMapGeneric_n = brighten(cMapGeneric,0.1);
    %b = bar(Pbar);
% don't change the color based on significance
%     for k=1:length(Pplot)
%         if Pplot(k)>0.05
%             cMapGeneric_n(k,:) = [235,235,235]/255; %brighten(cMapGeneric(k,:),0.95);
%         elseif Pplot(k)<=0.05 && Pplot(k)>=0.05/(length(Pplot))
%             cMapGeneric_n(k,:) = brighten(cMapGeneric(k,:),0.45);
%         end
%     end
    
%      ax{l}.XTick = 1:length(whatDiseases_Treatment_label);
%      ax{l}.XTickLabel = whatDiseases_Treatment_label(ix);
                
    %b.CData = cMapGeneric_n(ix,:);
    %b.FaceColor = 'flat';
    
    
    hold on;
    stem(Pbar, 'Marker','none', 'LineStyle',':', 'Color',[.25 .25 .25], 'LineWidth',2)
    
    b = scatter(1:length(ix), Pbar, 500, cMapGeneric_n(ix,:),...
        'MarkerFaceColor',[1 1 1],...
        'LineWidth',6);
    set(gcf, 'renderer', 'painters')
    xticks(1:length(ix));
    xticklabels(whatDiseases_Treatment_label(ix));
    xtickangle(90);
    xlim([0 5]); 
    ylim([0 4])
%     ax{l}.XTickLabelRotation = 90;
    xlabel('Disorder')
    ylabel('-log10(P)')    
    title(sprintf('%s',similarityTypes_label{l}),'interpreter','none')
        
        
    linkaxes([ax{:}],'y');
    set(gca,'FontSize', 18)
end




end