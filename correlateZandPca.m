doCoherence = 0;
doPlot = 0;
doRereference = 0;

zFolder = '/home/data/EEG/processed/Robi/zCorrelation';
pcaFolder = '/home/data/EEG/processed/Robi/coherencePca/';
zOutputFolder = '/media/eegDrive/zScoreTimecourse';
zPcaOutputFolder = '/media/eegDrive/zScorePca';
coherenceDataPath = '/media/eegDrive/';

files = dir(coherenceDataPath);
files([files.isdir]) = [];
filenameFilter = '_003';
if(length(filenameFilter) > 0)
    remove = cellfun(@length, strfind({files.name}, filenameFilter)) == 0;
end
files(remove) = [];

for fileCounter = 1:length(files)
    longFilename = files(fileCounter).name;
    shortFilename = longFilename(1:strfind(longFilename, '_63')-1);
    shortFilename = sprintf('%s.mat', shortFilename);
    fprintf('\n%s (%d of %d) %s', char(datetime), fileCounter, length(files), shortFilename);

    pcaBasePath = fullfile(pcaFolder, shortFilename);
    zPath = fullfile(zFolder, shortFilename);
    cohPath = fullfile(coherenceDataPath, longFilename); %'/media/eegDrive/ROBI_003_tx 3_630166729229512303coherenceStats.mat';
    zOutputPath = fullfile(zOutputFolder, shortFilename);
    zPcaPath = fullfile(zPcaOutputFolder, shortFilename);
    cohData = load(cohPath);
    badReferenceIndex = checkForBadReference(cohData.channelPairs(31).coherence);
    for i = 1:length(cohData.channelPairs)
        cohData.channelPairs(i).coherence(badReferenceIndex:end, :) = [];
        cohData.channelPairs(i).phaseAngle(badReferenceIndex:end, :) = [];
    end
    for i = 1:length(cohData.channels)
        cohData.channels(i).absolutePower(badReferenceIndex:end, :) = [];
        cohData.channels(i).relativePower(badReferenceIndex:end, :) = [];        
    end
    
    if(exist(zPath,'file'))
        zData = load(zPath);
    end
    if(exist(pcaBasePath,'file'))
        pcaBasePath = load(pcaBasePath);
    end
    zScoreTimecourse = zScore(cohData);
    save(zOutputPath, 'zScoreTimecourse');
    
    powColumns = size(zScoreTimecourse.absPower, 2) * size(zScoreTimecourse.absPower, 3);
    cohColumns = size(zScoreTimecourse.coherence, 2) * size(zScoreTimecourse.coherence, 3);
    timeRows = size(zScoreTimecourse.coherence, 1);
    timecourse = reshape(zScoreTimecourse.coherence, [timeRows, cohColumns]);
    timecourse(:,end+1:end+powColumns) = reshape(zScoreTimecourse.absPower, [timeRows, powColumns]);
    [dataPca.COEFF, dataPca.SCORE, dataPca.LATENT, dataPca.TSQUARED, dataPca.EXPLAINED, dataPca.MU] = pca(timecourse);
    maxCoefficient = 20;
    dataPca.COEFF(:, maxCoefficient+1:end) = [];
    dataPca.SCORE(:, maxCoefficient+1:end) = [];
    dataPca.LATENT(maxCoefficient+1:end) = [];
    dataPca.EXPLAINED(maxCoefficient+1:end) = [];
    dataPca.MU(maxCoefficient+1:end) = [];
    labelCounter = 1;
    freqLabels = {'delta','theta','alpha','beta','hibeta'};
    aniChanLocs = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2','F7','F8','T3','T4','T5','T6','Fz','Cz','Pz'};
    for i = 1:size(zScoreTimecourse.coherence,3)
        for j = 1:size(zScoreTimecourse.coherence,2)
            dataPca.pairLabels{labelCounter} = sprintf('%s %s', zScoreTimecourse.channelPairLabels{i}, freqLabels{j});
            labelCounter = labelCounter + 1;
        end
    end
    for i = 1:size(zScoreTimecourse.absPower,3)
        for j = 1:size(zScoreTimecourse.absPower,2)
            dataPca.pairLabels{labelCounter} = sprintf('%s %s abs', zScoreTimecourse.channelLabels{i}, freqLabels{j});
            labelCounter = labelCounter + 1;
        end
    end

    if(doPlot)
        plotCoherencePca(dataPca, dataPca.pairLabels, 1);
    end
    save(zPcaPath, 'dataPca');
end


if(false)
    zScoreFilename = '/home/data/EEG/data/ROBI/ROBI_003/tx 3/630166729229512303Events.txt';
    pcaFilename = '/media/eegDrive/zScorePca/ROBI_003_tx 3.mat';
    testScore = loadInterpolatedZScore(zScoreFilename);
    dataPca = load(pcaFilename);
    
    pcaData = dataPca.dataPca.SCORE(:,1);
    z = testScore.InstantZ';
    figure;
    hold on;
    plot(pcaData);
    plot(z, 'r');
    pan xon;
    zoom xon;
    
end

