function [coef_lasso, coef_linear, sigMeasures] = compare_lasso_coef(whatMeasures, whatDiseases_GWAS)
%X - predictor data: all gene scores from GWAS
%y - response: drug score
if nargin < 1
    whatMeasures = 'reduced'; % 'reduced' or 'all'; 
end

if nargin < 2
    whatDiseases_GWAS = {'ADHD', 'MDD2', 'SCZ', 'BIP2', 'DIABETES', 'HF'};
end

numGWAS = length(whatDiseases_GWAS);
coef_lasso = cell(numGWAS,4);
coef_linear = cell(numGWAS,4);
sigMeasures = cell(numGWAS,2);

% run lasso for each GWAS
for i = 1:numGWAS
    
    whatGWAS = whatDiseases_GWAS{i};
    [geneWeightsGWAS_ALL, drugScores_ord, measureNames] = give_GWASandDRUG_scores(whatGWAS, whatMeasures);
    
    % remove rows that are all NaN and measure names
    INDnan = find(all(isnan(geneWeightsGWAS_ALL),1));
    geneWeightsGWAS_ALL(:, INDnan) = [];
    measureNames(INDnan) = [];
    
    % apply linear regression
    mdl = fitlm(geneWeightsGWAS_ALL,drugScores_ord);
    coef_linear{i,1} = mdl.Coefficients(2:end,:); % don't take intercept
    coef_linear{i,2} = measureNames;
    A = anova(mdl,'summary'); 
    coef_linear{i,3} = A.F(2); 
    coef_linear{i,4} = A.pValue(2); 

    % find whihc p<0.05
    kk = find(mdl.Coefficients.pValue<0.05);
    if ~isempty(kk)
        sigMeasures{i,1} = measureNames(kk);
        sigMeasures{i,2} = mdl.Coefficients.pValue(kk);
    end
    
    % try lasso, in most cases, all values are set to 0, so not much use of that
    [B, FitInfo] = lasso(geneWeightsGWAS_ALL,drugScores_ord, 'CV',10, 'PredictorNames',measureNames);
    lassoPlot(B,FitInfo,'PlotType','CV');
    title(sprintf('%s',whatDiseases_GWAS{i}))
    legend('show')
    
    
    coef_lasso{i,1} = B(:,FitInfo.IndexMinMSE);
    coef_lasso{i,2} = B(:,50);
    coef_lasso{i,3} = measureNames;
    coef_lasso{i,4} = FitInfo.DF(FitInfo.IndexMinMSE);% number of parameters left at minSE 

    % try PLS
    geneWeightsGWAS_ALL(isnan(geneWeightsGWAS_ALL)) = 0;
    [XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(geneWeightsGWAS_ALL,drugScores_ord,size(geneWeightsGWAS_ALL,2));
    stats.W(:,1); 
    
    figure; plot(1:size(geneWeightsGWAS_ALL,2),cumsum(100*PCTVAR(2,:)),'-bo');
    xlabel('Number of PLS components');
    ylabel('Percent Variance Explained in y');
    title(sprintf('%s', whatGWAS))
end


%indPlot = find(tril(ones(numGWAS, numGWAS)));
whatRegressions = {'linear', 'lasso', 'PLS'};

for tt = 1:length(whatRegressions)
    whatRegression = whatRegressions{tt};
    r = nan(numGWAS, numGWAS);
    p = nan(numGWAS, numGWAS);
    %figure; set(gcf,'color','w');
    %i=1;
    for p1 = 1:numGWAS
        for p2 = p1:numGWAS
            
            %find measures that are common for each pair
            switch whatRegression
                case 'linear'
                    if p1~=p2
                        
                        [~, i1, i2] = intersect(coef_linear{p1,2}, coef_linear{p2,2});
                        [r(p1,p2), p(p1, p2)] = corr(coef_linear{p1,1}.Estimate(i1), coef_linear{p2,1}.Estimate(i2), 'rows', 'complete', 'type', 'Spearman');
                        
                    else
                        
                        r(p1,p2) = 0; % give a score of 0, so it doesn't disturb the scale; 
                        p(p1,p2) = coef_linear{p1,4};
                        
                    end
                    %subplot(numGWAS, numGWAS,indPlot(i));
                    %scatter(coef_linear{p1,1}.tStat(i1), coef_linear{p2,1}.tStat(i2));
                    %xlabel(whatDiseases_GWAS{p1})
                    %ylabel(whatDiseases_GWAS{p2})
                    %i=i+1;
                    
                case 'lasso'
                    
                    if p1~=p2
                        
                    [~, i1, i2] = intersect(coef_lasso{p1,3}, coef_lasso{p2,3});
                    [r(p1,p2), p(p1, p2)] = corr(coef_lasso{p1,2}(i1), coef_lasso{p2,2}(i2), 'rows', 'complete', 'type', 'Spearman');
%                     subplot(numGWAS, numGWAS,indPlot(i));
%                     scatter(coef_lasso{p1,2}(i1), coef_lasso{p2,2}(i2));
%                     xlabel(whatDiseases_GWAS{p1})
%                     ylabel(whatDiseases_GWAS{p2})
%                     i=i+1;
                    else
                        
                        r(p1,p2) = 0; % give a score of 0, so it doesn't disturb the scale; 
                        p(p1,p2) = coef_lasso{p1,4}; % this is a number of non-zero parameters at minSE
                        
                    end
                
                    
            end
        end
    end
    
    mVal = max(abs(r(:)));
    
    colors = cbrewer('div', 'RdBu', 64);
    colors = flipud(colors);
    figure; set(gcf,'color','w');
    imagesc(r); axis('square')
    colormap(colors);
    
    hold on;
    % plot p-values on top
    plot_matrixValues(p)
    
    yticks(1:numGWAS);
    yticklabels(whatDiseases_GWAS)
    
    xticks(1:numGWAS);
    xticklabels(whatDiseases_GWAS)
    
    caxis([-mVal mVal])
    colorbar
    title(sprintf('%s', whatRegression))
    
    
end

end
    
    
    