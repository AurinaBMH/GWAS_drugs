% run correlate_geneMeasures
clear all; close all; 

params = SetDefaultParams();
whatGWAS = params.whatGWAS; 
whatMeasures = 'allPsych'; 

for i=1:length(whatGWAS)
    
    correlate_geneMeasures(whatGWAS{i}, whatMeasures); 
    
    figureName = sprintf('figures/%s_geneMeasures', whatGWAS{i});
    print(gcf,figureName,'-dpng','-r300');
    
end
