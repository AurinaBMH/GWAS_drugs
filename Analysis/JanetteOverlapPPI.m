% For Janette
% A bit strange -> wants for
% - genes that overlap by annotation of SZP and BIP, and have protein interactions among them
% - genes that overlap by annotation of ADHD and ASD, and have protein interactions among them
% - genes that overlap by annotation of SZP and MDD, and, have protein interactions among them
% - genes that overlap by annotation of SZP, BIP and MDD, and, have protein interactions among them
%-------------------------------------------------------------------------------
% (i.e., get the list, then filter out genes that have no interactions in PPI network)
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
% Load data:
load(fullfile('Data','processedData.mat'),'Adj','eQTLproteinnames');
eQTLidentifier = importIdentifier();

%===============================================================================
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DOUBLES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%===============================================================================
% disease1 = 'SZP'; disease2 = 'BIP';
% disease1 = 'ADHD'; disease2 = 'ASD';
% disease1 = 'SZP'; disease2 = 'MDD';

isD1 = ismember(eQTLproteinnames,...
            unique(eQTLidentifier.Name(eQTLidentifier.(disease1) & ~eQTLidentifier.Partners)));
isD2 = ismember(eQTLproteinnames,...
            unique(eQTLidentifier.Name(eQTLidentifier.(disease2) & ~eQTLidentifier.Partners)));
isEither = (isD1|isD2);
numEither = sum(isEither);
fprintf(1,'%u genes have either (or both) %s/%s\n',numEither,disease1,disease2);
Adj_filter = Adj(isEither,isEither);

hasNoInteractions = (sum(Adj_filter)==0);
fprintf(1,'%u genes have no interactions within the PPI network -> %u\n',...
                sum(hasNoInteractions),numEither-sum(hasNoInteractions));

% Write out the original list:
writeOut(eQTLproteinnames,isEither,sprintf('%s_%s_either.csv',disease1,disease2));
% Now after filtering on P-P interactions:
writeOut(eQTLproteinnames(isEither),~hasNoInteractions,sprintf('%s_%s_either_PPIfilter.csv',disease1,disease2));

%===============================================================================
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TRIPLES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%===============================================================================
disease1 = 'SZP'; disease2 = 'BIP'; disease3 = 'MDD';

isD1 = ismember(eQTLproteinnames,...
            unique(eQTLidentifier.Name(eQTLidentifier.(disease1) & ~eQTLidentifier.Partners)));
isD2 = ismember(eQTLproteinnames,...
            unique(eQTLidentifier.Name(eQTLidentifier.(disease2) & ~eQTLidentifier.Partners)));
isD3 = ismember(eQTLproteinnames,...
            unique(eQTLidentifier.Name(eQTLidentifier.(disease3) & ~eQTLidentifier.Partners)));
isEither = (isD1|isD2|isD3);
numEither = sum(isEither);
fprintf(1,'%u genes annotated for at least one of: %s/%s/%s\n',numEither,disease1,disease2,disease3);
Adj_filter = Adj(isEither,isEither);

hasNoInteractions = (sum(Adj_filter)==0);
fprintf(1,'%u genes have no interactions within the PPI network -> %u\n',...
                sum(hasNoInteractions),numEither-sum(hasNoInteractions));

% Write out the original list:
writeOut(eQTLproteinnames,isEither,sprintf('%s_%s_%s_either.csv',disease1,disease2,disease3));
% Now after filtering on P-P interactions:
writeOut(eQTLproteinnames(isEither),~hasNoInteractions,sprintf('%s_%s_%s_either_PPIfilter.csv',disease1,disease2,disease3));
