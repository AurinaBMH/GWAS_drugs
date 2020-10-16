
whatDiseases_Treatment_ALL = params.whatDiseases_Treatment_ALL; 
[geneNamesDrug,drugScoresAll] = GiveMeNormalizedScoreVectors(whatDiseases_Treatment_ALL,'Drug');


figure('color', 'w'); 
for i=1:length(whatDiseases_Treatment_ALL)
    
    IND = find(drugScoresAll(:,i)); 
    N(i) = length(IND);
    V{i} = drugScoresAll(IND, i); 
    subplot(2,5,i); 
    histogram(V{i}); xlabel('Score'); ylabel('Number of targets'); 
    xlim([0 1])
    title(sprintf('%s, %d targets', whatDiseases_Treatment_ALL{i}, N(i)))
    G{i} = IND; 
    
end

% what is the overlap between gene targets across disorders? 
ND = length(whatDiseases_Treatment_ALL);
IN = zeros(ND); 
INperc = zeros(ND); 
for i=1:ND
    for j=1:ND
       
        IN(i,j) = length(intersect(G{i}, G{j})); 
        INperc(i,j) = length(intersect(G{i}, G{j}))/(length(unique(vertcat(G{i}, G{j})))); %(length(G{i})+length(G{j})); 
        
    end
end
[COL]=cbrewer('seq', 'Reds', 100);

IN(1:ND+1:end) = NaN; 
INperc(1:ND+1:end) = NaN;
figure('color', 'w'); 
subplot(1,2,1); imagesc(IN); axis 'square'; title('Number of shared targets'); 
colormap(COL); 
xticks(1:ND); xticklabels(whatDiseases_Treatment_ALL); 
yticks(1:ND); yticklabels(whatDiseases_Treatment_ALL); 
subplot(1,2,2); imagesc(INperc); axis 'square'; title('Proportion of shared targets')
colormap(COL); 
xticks(1:ND); xticklabels(whatDiseases_Treatment_ALL); 
yticks(1:ND); yticklabels(whatDiseases_Treatment_ALL); 
