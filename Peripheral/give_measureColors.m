function [colors, measureLabel, cluster] = give_measureColors(measureNames)

numMeasures = length(measureNames); 
cluster = zeros(length(measureNames),1); 
measureLabel = measureNames; 
 
colors = zeros(numMeasures, 3);
for tt=1:numMeasures
    if contains(measureNames{tt}, 'PPI') && contains(measureNames{tt}, 'num')
        colors(tt,:) = [31,120,180]; % dark blue
        A = split(measureNames{tt}, '.'); 
        measureLabel{tt} = A{1};
        cluster(tt) = 6; 
    elseif contains(measureNames{tt}, 'PPI') && contains(measureNames{tt}, 'perc')
        colors(tt,:) = [171,217,233]; % light blue
        A = split(measureNames{tt}, '.'); 
        measureLabel{tt} = A{1};
        cluster(tt) = 5; 
    elseif contains(measureNames{tt}, 'Allen')
        colors(tt,:) = [135,135,135]; % grey
        if contains(measureNames{tt}, 'mapped')
        measureLabel{tt} = 'AHBA_mapped'; 
        elseif contains(measureNames{tt}, 'eQTL')
        measureLabel{tt} = 'AHBA_eQTL'; 
        end
        cluster(tt) = 4; 
    elseif contains(measureNames{tt}, 'eQTL')
        colors(tt,:) = [178,223,138]; % light green
        A = split(measureNames{tt}, '.'); 
        measureLabel{tt} = A{1};
        cluster(tt) = 3; 
    elseif contains(measureNames{tt}, 'default')
        colors(tt,:) = [51,160,44]; % dark green
        A = split(measureNames{tt}, '.'); 
        measureLabel{tt} = A{1};
        cluster(tt) = 1; 
    elseif contains(measureNames{tt}, 'Combined')
        colors(tt,:) = [178,24,43];
        measureLabel{tt} = measureNames{tt}; 
        
    else % chromatin
        colors(tt,:) = [251,154,153]; % pink
        A = split(measureNames{tt}, '.'); 
        measureLabel{tt} = A{1};
        cluster(tt) = 2; 
        
    end
            
end

colors = (colors/255);

end