inputFolder = '/home/data/EEG/processed/Robi/coherencePca/rest';
inputFiles = dir(inputFolder);
inputFiles([inputFiles.isdir]) = [];

close all;

patientStarts = [1 3 7 11]
patientNumber = 3;

fileLimits = [patientStarts(patientNumber), patientStarts(patientNumber)+3];
startFile = fileLimits(1);
endFile = fileLimits(2);

fileIndices = startFile:endFile;
fileIndices = [8 10];

globalMin = realmax;
globalMax = realmin;
%for fileCounter = 1:length(inputFiles)
for fileCounter = fileIndices
    filename = fullfile(inputFolder, inputFiles(fileCounter).name);    
    a = load(filename);
    stdScores = std(a.cohPca.SCORE, 1);
    for i = 1:size(a.cohPca.COEFF, 2)
    normCoeff(:,i) = a.cohPca.COEFF(:,i) .* stdScores
    end
    
    data{fileCounter} = a;
    globalMin = min(globalMin, min(a.cohPca.COEFF(:,1)));
    globalMax = max(globalMax, max(a.cohPca.COEFF(:,1)));
end

for fileCounter = fileIndices  
    a = data{fileCounter};
    freqSize = size(a.cohPca.COEFF, 1) / 5;
    freqLabels = [{'delta'},{'theta'},{'alpha'},{'beta'},{'hibeta'}];
    
    for i = 1:5
        hasFreq = find(cellfun(@length, strfind(a.cohPca.labelChannelFreq, freqLabels{i})));
        labels = a.cohPca.labelChannelOnly(hasFreq);
        weights = a.cohPca.COEFF(hasFreq, 1);
        plotChannelPairs(labels, weights, [globalMin, globalMax], false);
        myTitle = sprintf('%s(%.1f)%s', freqLabels{i}, a.cohPca.EXPLAINED(1), inputFiles(fileCounter).name);
        myTitle = strrep(myTitle, '_', ' ');
        figureFolder = '/home/data/EEG/processed/Robi/figure';
        print(fullfile(figureFolder, sprintf('%s.png',myTitle)), '-dpng');        
        drawnow;        
    end    
end

tilefigs([5 length(fileIndices)]);