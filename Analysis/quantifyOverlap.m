function howMuchOverlap = quantifyOverlap(disease1,disease2,eQTLproteinnames,eQTLidentifier,Adj,doKeepPartners)

if nargin < 1
    disease1 = 'SZP';
end
if nargin < 2
    disease2 = 'BIP';
end
if nargin < 6
    doKeepPartners = false;
end

%-------------------------------------------------------------------------------

if doKeepPartners
    isD1 = ismember(eQTLproteinnames,...
                unique(eQTLidentifier.Name(eQTLidentifier.(disease1))));
    isD2 = ismember(eQTLproteinnames,...
                unique(eQTLidentifier.Name(eQTLidentifier.(disease2))));
else
    isD1 = ismember(eQTLproteinnames,...
                unique(eQTLidentifier.Name(eQTLidentifier.(disease1) & ~eQTLidentifier.Partners)));
    isD2 = ismember(eQTLproteinnames,...
                unique(eQTLidentifier.Name(eQTLidentifier.(disease2) & ~eQTLidentifier.Partners)));
end

isEither = (isD1|isD2);
fprintf(1,'%u proteins for %s or %s\n',sum(isEither),disease1,disease2);
Adj_filter = Adj(isEither,isEither);
isD1_filter = isD1(isEither);
isD2_filter = isD2(isEither);

offDiagonal = Adj_filter(isD1_filter,isD2_filter);
linkDensity = mean(offDiagonal(:));
howMuchOverlap = linkDensity;

end
