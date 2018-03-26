plotIndividual = 1;
filter = 'open';

files= getRobiDataFiles;

% inputFolder = '/media/eegDrive/zScoreTimecourse';
% files = dir(inputFolder);
% files([files.isdir]) = [];
filter = 'open';
files(find(cellfun(@length, strfind(files, filter)) == 0)) = [];

% inputFile = fullfile(inputFolder, 'ROBI_003_baseline eyes open.mat');
% exitFile = fullfile(inputFolder, 'ROBI_003_outcome eyes open.mat');
savedLocs = load('/home/data/EEG/processed/Robi/antChanlocs.mat');
allValues = zeros(1,length(files));
writeCounter = 1;
for fileCounter = 1:length(files)
    fprintf('\n%s: %d of %d', char(datetime), fileCounter, length(files));
    inputPath = files{fileCounter};
    slashes = strfind(inputPath, '/');
    outputFile = inputPath(slashes(end-2)+1:slashes(end)-1);
    outputFile = strrep(outputFile, '/','_');
    outputZFile = sprintf('/media/eegDrive/fast/fastZ/%s.mat', outputFile);
    outputFile = sprintf('/media/eegDrive/fast/%s.mat', outputFile);
    clear zScores;
    if(~exist(outputZFile,'file'))
        [inputFolder, inputFile, inputExt] = fileparts(inputPath);
        [hasNoise, noiseRatio] = checkForMainsNoise({[inputFile inputExt]}, inputFolder);
        if(~hasNoise)
            if(~exist(outputFile,'file'))
                data = loadRobiDataFile(files{fileCounter});
                EEG.data = (data(:, 1:33))';
                EEG.nbchan = 33;
                EEG.srate = 2048;
                EEG.chanlocs = savedLocs.chanlocs(1:33);
                [timeCourse.coherencePlot, timeCourse.times, timeCourse.powerPlot] = allChannelCoherence(EEG);
                timeCourse.basefile = files{fileCounter};
                timeCourse.mainsNoiseRatio = noiseRatio;
                timeCourse.badReferenceFramesDropped = length(timeCourse.times) - badReferenceCutoff;
                save(outputFile, 'timeCourse', '-v7.3');
            else
                load(outputFile);
            end
            [badReferenceCutoff] = checkForBadReference(timeCourse.coherencePlot(18).coherence(:,4));
            for i = 1:length(timeCourse.coherencePlot)
                timeCourse.coherencePlot(i).coherence(badReferenceCutoff:end,:) = [];
            end
            for i = 1:length(timeCourse.powerPlot)
                timeCourse.powerPlot(i).absolutePower(badReferenceCutoff:end,:) = [];
            end
        end
        zScores = zScore(timeCourse);
        save(outputZFile, 'zScores', '-v7.3');
    else
        if(exist(outputZFile, 'file'))
            load(outputZFile);
        end
    end
    if(exist('zScores','var'))
        meanZ.channelPairLabels = zScores.channelPairLabels;
        meanZ.absPower = mean(zScores.absPower,1);
        meanZ.relPower = mean(zScores.relPower,1);
        meanZ.coherence = mean(zScores.coherence,1);
        meanZ.basefile = files{fileCounter};
        meanZ.mainsNoiseRatio = noiseRatio;
        meanZ.badReferenceFramesDropped = length(timeCourse.times) - size(zScores.absPower, 1) - 1;% timeCourse.badReferenceFramesDropped;
        allZScores(writeCounter) = meanZ;
        fprintf('\nfileCounter: %d zScoreLength: %d',fileCounter, length(allZScores));
        writeCounter = writeCounter + 1;
    end
end
meanZScores.files = files;
meanZScores.zScores = allZScores;
save('/media/eegDrive/fast/fastZ/meanZScores.mat', 'allZScores', '-v7.3');