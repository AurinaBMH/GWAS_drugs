function [minLim,maxLim] = plot_nullDistribution(nullScores, rhos, plotWhere)

% nullScores - a vector of null scores that wre used to make the distribution
% rhos - a vector of matching values
% plotWhere - the place indicating where to plot the distribution: 
% -- when using pooled null, this should be the number of number of treatments+1,
% -- plotting on each bar, then this should be the index of the corresponding bar


[ff,x] = ksdensity(nullScores,linspace(min(nullScores),max(nullScores),500),'function','pdf');
ff = 0.4*ff/max(ff);

plot(ones(2,1)*(plotWhere)+ff,x,'k');
plot(ones(2,1)*(plotWhere)-ff,x,'k');


% Add horizontal lines to aid comparison to null:
null_p50 = quantile(nullScores,0.5);
%plot([1,plotWhere],ones(2,1)*null_p50,':k')
null_p95 = quantile(nullScores,0.95);
plot([plotWhere-0.5,plotWhere+0.5],ones(2,1)*null_p95,':k', 'LineWidth',2)
if range(rhos) > 0
    maxLim = max(rhos)*1.2;
    if maxLim < null_p95
        maxLim = null_p95*1.2;
    end
    minLim = min(rhos)*0.8;
    null_p10 = quantile(nullScores,0.1);
    if minLim > null_p10
        minLim = null_p10*0.8;
    end
    
end

end

