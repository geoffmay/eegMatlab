

function [  ] = testBatch3(  )
%TESTBATCH3 Summary of this function goes here
%   Detailed explanation goes here

% cd eeglab13_4_4b
% eeglab
% cd ..

inputFolder = '../data/';
%inputFolder = '/Users/Geoff/Box Sync/For Geoff Mays/PTSD MIND EEG Files';
outputFolder = '../processed/';

[tones, ant, other] = getBdfFiles();
% allFiles = [ant];
allFiles = cell(0);
files = dir(inputFolder);
for i = 1:length(files)
    file = files(i).name;
    if(strfind(file, '.bdf'))
        allFiles{length(allFiles)+1} = file;
    end
end

for i = length(allFiles):-1:1
    saveResults = strcat(outputFolder, allFiles{i});
    %if(~exist(strcat(saveResults, 'clusters.mat')))
        plotResults = false;
        smoothMaps = true;
        numberOfClusters = 400;
        bdfFileName = strcat(inputFolder, allFiles{i});
        [ topographicMaps, segmentCenters, segmentBoundaries, activationPlot, EEG ] = ...
            demoGFP(bdfFileName, numberOfClusters, smoothMaps, plotResults,saveResults );
    %end
end

end

