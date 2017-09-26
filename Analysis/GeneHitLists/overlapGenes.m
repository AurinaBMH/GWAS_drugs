function overlapGenes(disease1,disease2)

if nargin < 1
    disease1 = 'SZP';
end
if nargin < 2
    disease2 = 'BIP';
end

%-------------------------------------------------------------------------------
% Load data:
load(fullfile('Data','processedData_eQTL.mat'),'Adj','proteinNames');
eQTLidentifier = importIdentifier();

%-------------------------------------------------------------------------------

isD1 = ismember(proteinNames,...
            unique(eQTLidentifier.Name(eQTLidentifier.(disease1) & ~eQTLidentifier.Partners)));
isD2 = ismember(proteinNames,...
            unique(eQTLidentifier.Name(eQTLidentifier.(disease2) & ~eQTLidentifier.Partners)));

isBoth = (isD1&isD2);

% Write full eQTL lists for each:
writeOut(proteinNames,isD1,sprintf('%s_eQTL.csv',disease1));
writeOut(proteinNames,isD2,sprintf('%s_eQTL.csv',disease2));

%-------------------------------------------------------------------------------
% List proteins overlapping:
%-------------------------------------------------------------------------------
fprintf(1,'%u overlapping annotations:\n',sum(isBoth));
writeOut(proteinNames,isBoth,sprintf('%s_%s_overlapping.csv',disease1,disease2));
overLappingNames = proteinNames(isBoth);
disp(overLappingNames);

%-------------------------------------------------------------------------------
% Filter all to the subset of genes involved in either disorder
%-------------------------------------------------------------------------------
isEither = (isD1|isD2);
numEither = sum(isEither);
isBoth_filter = isBoth(isEither);
isD1_filter = isD1(isEither);
isD1spec_filter = isD1_filter & ~isBoth_filter;
isD2_filter = isD2(isEither);
isD2spec_filter = isD2_filter & ~isBoth_filter;
Adj_filter = Adj(isEither,isEither);
proteinNames_filter = proteinNames(isEither);

%-------------------------------------------------------------------------------
% Remove self-interactions
%-------------------------------------------------------------------------------
fprintf(1,'EXCLUDING %u self-interactions\n',sum(diag(Adj_filter)));
Adj_filter(logical(eye(size(Adj_filter)))) = 0;

%-------------------------------------------------------------------------------
% crossFilt tells us where D1 interacts with D2
%-------------------------------------------------------------------------------
doesCrossInteract = countUnique(isD1_filter,isD2_filter,Adj_filter);
fprintf(1,'%u Interacting protein-coding genes:\n',sum(doesCrossInteract));
writeOut(proteinNames_filter,doesCrossInteract,...
                    sprintf('%s_%s_interact_any.csv',disease1,disease2));
% List proteins overlapping:
% interactingNames = proteinNames_filter(doesCrossInteract);
% disp(interactingNames);

%-------------------------------------------------------------------------------
% List proteins overlapping and disease1; disease2:
%-------------------------------------------------------------------------------
crossInteract_d1 = doesCrossInteract & isD1_filter;
crossInteract_d2 = doesCrossInteract & isD2_filter;

writeOut(proteinNames_filter,crossInteract_d1,...
            sprintf('%s_%s_interact_is%s.csv',disease1,disease2,disease1));
writeOut(proteinNames_filter,crossInteract_d2,...
            sprintf('%s_%s_interact_is%s.csv',disease1,disease2,disease2));


fprintf(1,'%u %s-annotated and interacting protein-coding genes:\n',...
                        sum(crossInteract_d1),disease1);
% theNames = proteinNames(fisEither(crossInteract_d1));
fprintf(1,'%u %s-annotated and interacting protein-coding genes:\n',...
                        sum(crossInteract_d2),disease2);
% theNames = proteinNames(fisEither(crossInteract_d2));
% for i = 1:length(theNames)
%     fprintf(1,'%s\n',theNames{i});
% end

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


%-------------------------------------------------------------------------------
% Write:
writeOut(proteinNames_filter,uniqueD1,...
            sprintf('%s_%s_unique_%s.csv',disease1,disease2,disease1));
writeOut(proteinNames_filter,uniqueD2,...
            sprintf('%s_%s_unique_%s.csv',disease1,disease2,disease2));

end
