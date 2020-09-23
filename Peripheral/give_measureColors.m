function colors = give_measureColors(measureNames)

numMeasures = length(measureNames); 

colors = zeros(numMeasures, 3);
for tt=1:numMeasures
    if contains(measureNames{tt}, 'PPI') && contains(measureNames{tt}, 'num')
        colors(tt,:) = [31,120,180]; % dark blue
    elseif contains(measureNames{tt}, 'PPI') && contains(measureNames{tt}, 'perc')
        colors(tt,:) = [171,217,233]; % light blue
    elseif contains(measureNames{tt}, 'Allen')
        colors(tt,:) = [135,135,135]; % grey
    elseif contains(measureNames{tt}, 'eQTL')
        colors(tt,:) = [178,223,138]; % light green
    elseif contains(measureNames{tt}, 'default')
        colors(tt,:) = [51,160,44]; % dark green
    elseif contains(measureNames{tt}, 'Combined')
        colors(tt,:) = [178,24,43];
    else % chromatin
        colors(tt,:) = [251,154,153]; % pink
    end
end

colors = (colors/255);

end