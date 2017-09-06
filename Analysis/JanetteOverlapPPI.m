% For Janette
% A bit strange -> wants for
% - genes that overlap by annotation of SZP and BIP, and have protein interactions among them
% - genes that overlap by annotation of ADHD and ASD, and have protein interactions among them
% - genes that overlap by annotation of SZP and MDD, and, have protein interactions among them
% - genes that overlap by annotation of SZP, BIP and MDD, and, have protein interactions among them
%-------------------------------------------------------------------------------
% (i.e., get the list, then filter out genes that have no interactions in PPI network)
%-------------------------------------------------------------------------------
doeQTL = true;

%-------------------------------------------------------------------------------
% Load data:
%-------------------------------------------------------------------------------
if doeQTL
    dataFile = 'processedData_eQTL.mat';
    geneIdentifier = importIdentifier();
else
    dataFile = 'processedData_Mapped.mat';
    geneIdentifier = ImportIdentifierMapped();
end
load(fullfile('Data','processedData_eQTL.mat'),'Adj','proteinNames');

%===============================================================================
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DOUBLES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%===============================================================================
% disease1 = 'SZP'; disease2 = 'BIP';
% disease1 = 'ADHD'; disease2 = 'ASD';
% disease1 = 'SZP'; disease2 = 'MDD';

isD1 = ismember(proteinNames,...
            unique(geneIdentifier.Name(geneIdentifier.(disease1) & ~geneIdentifier.Partners)));
isD2 = ismember(proteinNames,...
            unique(geneIdentifier.Name(geneIdentifier.(disease2) & ~geneIdentifier.Partners)));
isEither = (isD1|isD2);
numEither = sum(isEither);
fprintf(1,'%u genes have either (or both) %s/%s\n',numEither,disease1,disease2);
Adj_filter = Adj(isEither,isEither);

hasNoInteractions = (sum(Adj_filter)==0);
fprintf(1,'%u genes have no interactions within the PPI network -> %u\n',...
                sum(hasNoInteractions),numEither-sum(hasNoInteractions));

% Write out the original list:
writeOut(proteinNames,isEither,sprintf('%s_%s_either.csv',disease1,disease2));
% Now after filtering on P-P interactions:
writeOut(proteinNames(isEither),~hasNoInteractions,sprintf('%s_%s_either_PPIfilter.csv',disease1,disease2));

%===============================================================================
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TRIPLES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%===============================================================================
disease1 = 'SZP'; disease2 = 'BIP'; disease3 = 'MDD';

isD1 = ismember(proteinNames,...
            unique(geneIdentifier.Name(geneIdentifier.(disease1) & ~geneIdentifier.Partners)));
isD2 = ismember(proteinNames,...
            unique(geneIdentifier.Name(geneIdentifier.(disease2) & ~geneIdentifier.Partners)));
isD3 = ismember(proteinNames,...
            unique(geneIdentifier.Name(geneIdentifier.(disease3) & ~geneIdentifier.Partners)));
isEither = (isD1|isD2|isD3);
numEither = sum(isEither);
fprintf(1,'%u genes annotated for at least one of: %s/%s/%s\n',numEither,disease1,disease2,disease3);
Adj_filter = Adj(isEither,isEither);

hasNoInteractions = (sum(Adj_filter)==0);
fprintf(1,'%u genes have no interactions within the PPI network -> %u\n',...
                sum(hasNoInteractions),numEither-sum(hasNoInteractions));

% Write out the original list:
writeOut(proteinNames,isEither,sprintf('%s_%s_%s_either.csv',disease1,disease2,disease3));
% Now after filtering on P-P interactions:
writeOut(proteinNames(isEither),~hasNoInteractions,sprintf('%s_%s_%s_either_PPIfilter.csv',disease1,disease2,disease3));
