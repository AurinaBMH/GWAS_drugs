function drugClassTable = DrugClassImport()
% Idea is to assign each drug to a class
%-------------------------------------------------------------------------------
fid = fopen('8_2_drug_ATC_matrix.csv','r');
C = textscan(fid,'%s%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u','Delimiter',',','HeaderLines',1);
fclose(fid);

drugName = C{1};
isA = C{2};
isB = C{3};
isC = C{4};
isD = C{5};
isG = C{6};
isH = C{7};
isJ = C{8};
isL = C{9};
isM = C{10};
isN = C{11};
isP = C{12};
isR = C{13};
isS = C{14};
isV = C{15};
isNil = C{16};

classes = {'A','B','C','D','G','H','J','L','M','N','P','R','S','V','Nil'};
isClass = logical([horzcat(C{2:16})]);
numDrugs = length(drugName);

whatClass = cell(numDrugs,1);
for i = 1:numDrugs
    whatClass{i} = classes{isClass(i,:)};
end
whatClass = categorical(whatClass);

drugClassTable = table(drugName,whatClass);

fprintf(1,'Classifications for %u drugs\n',height(drugClassTable));

end
