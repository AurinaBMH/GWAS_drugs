function [diseaseResultsR, diseaseResultsP, similarityTypes,ALLmeasures] = compareGWASvsDRUGmatches(whatDiseases_GWAS, whatNull, Dname, PPImeasures_names, similarityTypes)

if nargin <2
    whatNull = 'randomDisease';
end

if nargin <3
    % if the target drug list not provided, compare with itself
    Dname = whatDiseases_GWAS{1};
end

if nargin <4
    PPImeasures_names = {'numPPIneighbors1','percPPIneighbors1','gwasPPIneighbors1', 'weiPPIneighbors1','expWeiPPIneighbors1', 'weigwasPPIneighbors1', 'numPPIneighbors2','percPPIneighbors2', 'gwasPPIneighbors2', 'weiPPIneighbors2','expWeiPPIneighbors2', 'weigwasPPIneighbors2'};
end

if nargin<5
    similarityTypes = {'Adult_brain', 'AllenMeanCoexpMapped', 'AllenMeanCoexpeQTLbrain', ...
        'Astro', 'Fetal_brain', 'MAGMAdefault', 'Neuro', ...
        'PPI_eQTLbrain_th0', 'PPI_eQTLbrain_th400', 'PPI_eQTLbrain_th600', 'PPI_eQTLbrain_th900', ...
        'PPI_mapped_th0', 'PPI_mapped_th400', 'PPI_mapped_th600', 'PPI_mapped_th900', ...
        'eQTLHeart_Left_Ventricle', 'eQTLLiver', 'eQTLWhole_Blood', 'eQTLbrain'};
end

if strcmp(Dname, 'HF')
    Dname = 'cardiology';
end
Dname = Dname(isstrprop(Dname,'alpha'));

% measures for non PPI mappings
OTHERmeasures_names = {'P', 'r'};
whatThreshold = 'BF';

PPInum = length(PPImeasures_names);
OTHERnum = length(OTHERmeasures_names);

diseaseResultsR = nan(length(similarityTypes), PPInum+OTHERnum);
diseaseResultsP = nan(length(similarityTypes), PPInum+OTHERnum);
diseaseResultsISbest = strings(length(similarityTypes), PPInum+OTHERnum);

for s=1:length(similarityTypes)
    
    if contains(similarityTypes{s},'PPI')
        for p=1:PPInum
            
            whatProperty = PPImeasures_names{p};
            [rhos ,pVals, whatDiseases_Treatment] = DistinguishingCharBar(similarityTypes{s}, whatProperty, whatNull, whatThreshold, whatDiseases_GWAS, false);
            % select rho and p values for a selected disorder
            
            takeVal = contains(whatDiseases_Treatment, Dname, 'IgnoreCase',true);
            
            diseaseResultsR(s,p+OTHERnum) = rhos(takeVal);
            diseaseResultsP(s,p+OTHERnum) = pVals(takeVal);
            
            % print on screen if NaN where shouldn't and replace it's
            % p-value with 1;
            if isnan(rhos(takeVal))
                warning('Dot product is NaN, replacing p-val with 1\n')
                diseaseResultsP(s,p+OTHERnum) = 1;
            end
            % test if the selected p-val is the most significant
            if diseaseResultsP(s,p+OTHERnum)<0.05 && diseaseResultsP(s,p+OTHERnum)<=min(pVals)
            %if diseaseResultsP(s,p+OTHERnum)<0.05 && diseaseResultsR(s,p+OTHERnum)>=max(rhos)
                diseaseResultsISbest(s,p+OTHERnum) = 'x';
            end
        end
    else % for non-PPI measures
        
        if ~contains(similarityTypes{s},'Allen')
            whatProperty = 'P';
        elseif contains(similarityTypes{s},'Allen')
            whatProperty = 'r';
        end
        
        [rhos ,pVals, whatDiseases_Treatment] = DistinguishingCharBar(similarityTypes{s}, whatProperty, whatNull, whatThreshold, whatDiseases_GWAS, false);
        takeVal = contains(whatDiseases_Treatment, Dname, 'IgnoreCase',true);
        
        colN = find(strcmp(OTHERmeasures_names,whatProperty));
        
        diseaseResultsR(s,colN) = rhos(takeVal);
        diseaseResultsP(s,colN) = pVals(takeVal);
        % significant and smallest p-value overall
        if pVals(takeVal)<0.05 && pVals(takeVal)<=min(rhos)
        %if pVals(takeVal)<0.05 && rhos(takeVal)>=max(rhos)
            diseaseResultsISbest(s,colN) = 'x';
        end
    end
    
end

ALLmeasures = horzcat(OTHERmeasures_names, PPImeasures_names);
similarityTypes_label = strrep(similarityTypes(:),'_',' '); % remove _ for plotting
colors = cbrewer('seq', 'Reds', 64);

figure; set(gcf,'color','w');
set(gcf, 'Position', [500 500 750 750])

imagesc(diseaseResultsR); axis('square')

%colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
hold on;
% plot p-values on top
plot_matrixValues(diseaseResultsP)
hold on;
plot_matrixValues(diseaseResultsISbest)
hold off


yticks(1:length(similarityTypes));
yticklabels(similarityTypes_label)

xticks(1:PPInum+OTHERnum);
xticklabels(ALLmeasures)

xtickangle(45);
box off

title(sprintf('%s GWAS vs %s Drugs \n null=%s', whatDiseases_GWAS{1}, Dname, whatNull))
colorbar
caxis([0 0.4])

end
