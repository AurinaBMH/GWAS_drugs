function howMuchOverlap = quantifyOverlap(disease1,disease2,eQTLproteinnames,eQTLidentifier,Adj)

if nargin < 1
    disease1 = 'SZP';
end
if nargin < 2
    disease2 = 'BIP';
end

%-------------------------------------------------------------------------------

isD1 = ismember(eQTLproteinnames,...
            unique(eQTLidentifier.Name(eQTLidentifier.(disease1) & ~eQTLidentifier.Partners)));
isD2 = ismember(eQTLproteinnames,...
            unique(eQTLidentifier.Name(eQTLidentifier.(disease2) & ~eQTLidentifier.Partners)));

isEither = (isD1|isD2);
Adj_filter = Adj(isEither,isEither);
isD1_filter = isD1(isEither);
isD2_filter = isD2(isEither);

offDiagonal = Adj_filter(isD1_filter,isD2_filter);
linkDensity = mean(offDiagonal(:));
howMuchOverlap = linkDensity;

end
