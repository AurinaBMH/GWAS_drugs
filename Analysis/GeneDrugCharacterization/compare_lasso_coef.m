function [coef_lasso, coef_linear, sigMeasures] = compare_lasso_coef(whatDiseases_GWAS)
%X - predictor data: all gene scores from GWAS
%y - response: drug score
if nargin < 1
    whatDiseases_GWAS = {'ADHD', 'MDD2', 'SCZ', 'BIP2', 'DIABETES', 'HF'};
end

numGWAS = length(whatDiseases_GWAS);
coef_lasso = cell(numGWAS,3);
coef_linear = cell(numGWAS,2);
sigMeasures = cell(numGWAS,2);

% run lasso for each GWAS
for i = 1:numGWAS
    
    whatGWAS = whatDiseases_GWAS{i};
    [geneWeightsGWAS_ALL, drugScores_ord, measureNames] = give_GWASandDRUG_scores(whatGWAS);
    
    % remove rows that are all NaN and measure names
    INDnan = find(all(isnan(geneWeightsGWAS_ALL),1));
    geneWeightsGWAS_ALL(:, INDnan) = [];
    measureNames(INDnan) = [];
    
    % apply linear regression
    mdl = fitlm(geneWeightsGWAS_ALL,drugScores_ord);
    coef_linear{i,1} = mdl.Coefficients(2:end,:); % don't take intercept
    coef_linear{i,2} = measureNames;
    % find whihc p<0.05
    kk = find(mdl.Coefficients.pValue<0.05);
    if ~isempty(kk)
        sigMeasures{i,1} = measureNames(kk);
        sigMeasures{i,2} = mdl.Coefficients.pValue(kk);
    end
    
    % try lasso, in most cases, all values are set to 0, so not much use of that
    [B, FitInfo] = lasso(geneWeightsGWAS_ALL,drugScores_ord, 'CV',10);
    lassoPlot(B,FitInfo,'PlotType','CV');
    idxLambda1SE = FitInfo.Index1SE;
    idxLambdaMinMSE = FitInfo.IndexMinMSE;
    coef_lasso{i,1} = B(:,idxLambda1SE);
    coef_lasso{i,2} = B(:,idxLambdaMinMSE);
    coef_lasso{i,3} = measureNames;
    
    legend('show') % Show legend
end


indPlot = find(tril(ones(numGWAS, numGWAS),-1));
whatRegressions = {'linear', 'lasso'};

for tt = 1:length(whatRegressions)
    whatRegression = whatRegressions{tt};
    r = nan(numGWAS, numGWAS);
    figure; set(gcf,'color','w');
    i=1;
    for p1 = 1:numGWAS
        for p2 = p1+1:numGWAS
            
            %find measures that are common for each pair
            switch whatRegression
                case 'linear'
                    [~, i1, i2] = intersect(coef_linear{p1,2}, coef_linear{p2,2});
                    r(p1,p2) = corr(coef_linear{p1,1}.tStat(i1), coef_linear{p2,1}.tStat(i2), 'rows', 'complete', 'type', 'Spearman');
                    subplot(numGWAS, numGWAS,indPlot(i));
                    scatter(coef_linear{p1,1}.tStat(i1), coef_linear{p2,1}.tStat(i2));
                    xlabel(whatDiseases_GWAS{p1})
                    ylabel(whatDiseases_GWAS{p2})
                    i=i+1;
                    
                case 'lasso'
                    
                    [~, i1, i2] = intersect(coef_lasso{p1,3}, coef_lasso{p2,3});
                    r(p1,p2) = corr(coef_lasso{p1,2}(i1), coef_lasso{p2,2}(i2), 'rows', 'complete', 'type', 'Spearman');
                    subplot(numGWAS, numGWAS,indPlot(i));
                    scatter(coef_lasso{p1,2}(i1), coef_lasso{p2,2}(i2));
                    xlabel(whatDiseases_GWAS{p1})
                    ylabel(whatDiseases_GWAS{p2})
                    i=i+1;
                    
            end
        end
    end
    
    mVal = max(abs(r(:)));
    
    colors = cbrewer('div', 'RdBu', 64);
    colors = flipud(colors);
    figure; set(gcf,'color','w');
    imagesc(r); axis('square')
    colormap(colors);
    
    yticks(1:numGWAS);
    yticklabels(whatDiseases_GWAS)
    
    xticks(1:numGWAS);
    xticklabels(whatDiseases_GWAS)
    
    caxis([-mVal mVal])
    colorbar
    title(sprintf('%s', whatRegression))
    
    
end

end
    
    
    