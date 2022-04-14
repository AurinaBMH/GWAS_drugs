% GenerateResultsTables
% Idea is to run through pipeline for each disease, and save results tables
% for each:
function GenerateResultsTables(whichYear)

if nargin < 1
    whichYear = '2022';
end

% Diseases with SNP info, etc.:
params = SetDefaultParams();
if strcmp(whichYear, '2022')
    whatDiseases = params.whatGWAS_2022;
elseif strcmp(whichYear, '2021')
    whatDiseases = params.whatGWAS_2021;
end

whatThreshold='BF';

% load information on genes mapped using different methods
load(sprintf('DataOutput_2022/GWAS_disordersMAGMA_%s.mat', whichYear))

%-------------------------------------------------------------------------------
numDiseases = length(whatDiseases);
for k = 1:numDiseases
    %for k = 1:numDiseases
    whatDisease = whatDiseases{k};
    geneScores = pipeline(DISORDERlist, whatDisease, whatThreshold);
    % Save:
    fileName = sprintf('resultsTable_%s_%s_%s_%s_drugbank.mat',whatDisease, whatThreshold, params.whatDrugTargets, params.whatTargets);
    fileName = fullfile('DataOutput_2022',fileName);
    save(fileName,'geneScores');
    fprintf(1,'Saved results to %s!!!\n\n\n',fileName);
end
end
