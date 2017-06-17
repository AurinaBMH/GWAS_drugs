
disease1 = 'SZP';
disease2 = 'BIP';

isD1 = ismember(eQTLproteinnames,...
            unique(eQTLidentifier.Name(eQTLidentifier.(disease1) & ~eQTLidentifier.Partners)));
isD2 = ismember(eQTLproteinnames,...
            unique(eQTLidentifier.Name(eQTLidentifier.(disease2) & ~eQTLidentifier.Partners)));

isBoth = (isD1&isD2);
% List proteins overlapping:
fprintf(1,'%u overlapping annotations:\n',sum(isBoth));
overLappingNames = eQTLproteinnames(isBoth);
disp(overLappingNames);

f_takeUpper = @(x) triu(x);

% Interactions
isEither = (isD1|isD2);
fisEither = find(isEither);
numEither = sum(isEither);
isBoth_filter = isBoth(isEither);
isD1_filter = isD1(isEither);
isD1spec_filter = isD1_filter & ~isBoth_filter;
isD2_filter = isD2(isEither);
isD2spec_filter = isD2_filter & ~isBoth_filter;
Adj_filter = Adj(isEither,isEither);
% crossFilt tells us where D1 interacts with D2
crossFilt = repmat(isD1_filter,1,numEither) & repmat(isD2_filter,1,numEither)';
crossFilt = (crossFilt | crossFilt');
crossInteraction = (Adj_filter & crossFilt);
doesCrossInteract = any(crossInteraction,1);
interactingNames = eQTLproteinnames(fisEither(doesCrossInteract));
% List proteins overlapping:
fprintf(1,'%u Interacting protein-coding genes:\n',length(interactingNames));
disp(interactingNames);

% List proteins overlapping and BIP:
crossInteract_d2 = any(crossInteraction,2) & isD2_filter;
fprintf(1,'%u %s-annotated and interacting protein-coding genes:\n',...
                        sum(crossInteract_d2),disease2);
theNames = eQTLproteinnames(fisEither(crossInteract_d2));
for i = 1:length(theNames)
    fprintf(1,'%s\n',theNames{i});
end

crossInteract_d1 = any(crossInteraction,2) & isD1_filter;
fprintf(1,'%u %s-annotated and interacting protein-coding genes:\n',...
                        sum(crossInteract_d1),disease1);
theNames = eQTLproteinnames(fisEither(crossInteract_d1));
for i = 1:length(theNames)
    fprintf(1,'%s\n',theNames{i});
end

%-------------------------------------------------------------------------------
% Genes involved in specific types of cross-disorder interactions
%-------------------------------------------------------------------------------
% D1/D2 <-> D1
doesBothD1 = countUnique(isBoth_filter,isD1spec_filter,Adj_filter);
fprintf(1,'%u %s/%s - %s\n',sum(doesBothD1),disease1,disease2,disease1);

% D1/D2 <-> D2
doesBothD2 = countUnique(isBoth_filter,isD2spec_filter,Adj_filter);
fprintf(1,'%u %s/%s - %s\n',sum(doesBothD2),disease1,disease2,disease2);

% D1/D2 <-> D1/D2
doesBothBoth = countUnique(isBoth_filter,isBoth_filter,Adj_filter);
fprintf(1,'%u %s/%s - %s/%s\n',sum(doesBothBoth),disease1,disease2,disease1,disease2);

% D1 <-> D2
doesD1D2 = countUnique(isD1spec_filter,isD2spec_filter,Adj_filter);
fprintf(1,'%u %s - %s\n',sum(doesD1D2),disease1,disease2);

% D1 <-> []
uniqueD1 = countUnique(isD1spec_filter,~isD2_filter,Adj_filter);
fprintf(1,'%u %s - []\n',sum(uniqueD1),disease1);

% D2 <-> []
uniqueD2 = countUnique(isD2spec_filter,~isD1_filter,Adj_filter);
fprintf(1,'%u %s - <>\n',sum(uniqueD2),disease2);

d2overlapList = (doesD1D2 | doesBothBoth | doesBothD2) & isD2_filter;
