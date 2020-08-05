function compare_lasso_coef(whatDiseases_GWAS, whatNull)
%X - predictor data: all gene scores from GWAS
%y - response: drug score

[geneNamesDrug,drugScoresAll] = GiveMeNormalizedScoreVectors(whatDiseases_Treatment,'Drug');


% creat a vector of rho values
B = lasso(X,y) 
end