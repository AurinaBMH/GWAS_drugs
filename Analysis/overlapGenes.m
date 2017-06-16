
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


both_D1 = repmat(isBoth_filter,1,numEither) & repmat(isD1spec_filter,1,numEither)';
both_D1 = (both_D1|both_D1');
cross_both_d1 = (Adj_filter & both_D1);
doesBothD1 = any(cross_both_d1,2);
fprintf(1,'%u %s/%s - %s\n',sum(doesBothD1),disease1,disease2,disease1);

both_D2 = repmat(isBoth_filter,1,numEither) & repmat(isD2spec_filter,1,numEither)';
both_D2 = (both_D2|both_D2');
cross_both_d2 = (Adj_filter & both_D2);
doesBothD2 = any(cross_both_d2,1);
fprintf(1,'%u %s/%s - %s\n',sum(doesBothD2),disease1,disease2,disease2);

both_both = repmat(isBoth_filter,1,numEither) & repmat(isBoth_filter,1,numEither)';
both_both = (both_both|both_both');
cross_both_both = (Adj_filter & both_both);
doesBothBoth = any(cross_both_both,2);
fprintf(1,'%u %s/%s - %s/%s\n',sum(doesBothBoth),disease1,disease2,disease1,disease2);

D1_D2 = repmat(isD1spec_filter,1,numEither) & repmat(isD2spec_filter,1,numEither)';
D1_D2 = (D1_D2|D1_D2');
cross_d1_d2 = (Adj_filter & D1_D2);
doesD1D2 = any(cross_d1_d2,2);
fprintf(1,'%u %s - %s\n',sum(doesD1D2),disease1,disease2);

d2overlapList = (doesD1D2 | doesBothBoth | doesBothD2) & isD2_filter;
