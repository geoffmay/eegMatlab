%plot and save maps for each ica component

clear;

percentileThreshold = 10;

inputFolder = '/home/data/EEG/processed/Robi/coherenceIca/old2/';
outputFolder = '/media/eegDrive/figs/robi003BaselineMastoid ICA components/';

inputFile = 'ROBI_003_baseline eyes open-MastoidsFast.mat';
load(fullfile(inputFolder, inputFile));
%load('/home/data/EEG/processed/Robi/coherenceIca/old2/ROBI_003_baseline eyes open-MastoidsFast.mat')
a = reshape(cohIca.separatingMatrix, [1 numel(cohIca.separatingMatrix)]);
b = sort(a);
coherenceIndexRange = 1:(30*29/2*5);
for i = 1:size(cohIca.separatingMatrix, 1)
    values = cohIca.separatingMatrix(i,:);
    labels = cohIca.labelChannelFreq;
    lowerThreshold = prctile(values(coherenceIndexRange), percentileThreshold);
    upperThreshold = prctile(values(coherenceIndexRange), 100-percentileThreshold);
    withinThreshold = values(coherenceIndexRange) > lowerThreshold & values(coherenceIndexRange) < upperThreshold;
    labels(withinThreshold) = [];
    values(withinThreshold) = [];
    
    outputFile = strrep(inputFile, '.mat', '');
    outputFile = sprintf('%sComponent%04d.png',outputFile,i);
    outputPath = fullfile(outputFolder, outputFile);    
    close all;
    fig = plotCoherencePca(values, labels);
    print(outputPath, '-dpng', '-r0');
end
