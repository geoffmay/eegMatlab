doCoherence = 0;
doPlot = 0;
doRereference = 0;

% zFolder = '/home/data/EEG/processed/Robi/zCorrelation';
% pcaFolder = '/home/data/EEG/processed/Robi/coherencePca/';
zOutputFolder = '/media/eegDrive/zScoreTimecourse';
% zPcaOutputFolder = '/media/eegDrive/zScorePca';
coherenceDataPath = '/media/eegDrive/';

files = dir(coherenceDataPath);
files([files.isdir]) = [];
filenameFilter = 'closed';
if(length(filenameFilter) > 0)
    remove = cellfun(@length, strfind({files.name}, filenameFilter)) == 0;
end
files(remove) = [];

for fileCounter = 1:length(files)
    longFilename = files(fileCounter).name;
    shortFilename = longFilename(1:strfind(longFilename, '_63')-1);
    shortFilename = sprintf('%s.mat', shortFilename);
    fprintf('\n%s (%d of %d) %s', char(datetime), fileCounter, length(files), shortFilename);
    
    %     pcaBasePath = fullfile(pcaFolder, shortFilename);
    %     zPath = fullfile(zFolder, shortFilename);
    cohPath = fullfile(coherenceDataPath, longFilename); %'/media/eegDrive/ROBI_003_tx 3_630166729229512303coherenceStats.mat';
    zOutputPath = fullfile(zOutputFolder, shortFilename);
    %     zPcaPath = fullfile(zPcaOutputFolder, shortFilename);
    if(~exist(zOutputPath, 'file'))
        cohData = load(cohPath);
        zScoreTimecourse = zScore(cohData);
        save(zOutputPath, 'zScoreTimecourse');
    end
end