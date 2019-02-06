folder = 'C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\phaseSlope';
outputFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\convolvedPhaseSlope';
eegFramesPerSecond = 1000;
phaseSlopeFramesPerSecond = 10;

mriFiles = fileLocations;
mriFiles = mriFiles.ghermanFiles;
eegFilenames = dir(folder);
eegFilenames([eegFilenames.isdir]) = [];
if(~exist(outputFolder, 'file'))
    mkdir(outputFolder);
end
for i = 1:length(eegFilenames)
    filename = eegFilenames(i).name;
    filename = strrep(filename, 'EEG_data_sub-', '');
    subjectNumber = filename(1:2);
    runNumber = filename(8:9);
    
    mriFolder = 'C:\Users\Neuro\Documents\MATLAB\data\GhermanPhilastides\EEG_events_volumemarkers\new_EEG_events_volumemarkers\';
    mriFilename = fullfile(mriFolder, sprintf('EEG_events_sub-%s_run-%s.mat', subjectNumber, runNumber));
    mriFilenames{i} = mriFilename;
end

fileCounter = 1;

while(fileCounter <= length(eegFilenames))
    filename = fullfile(folder, eegFilenames(fileCounter).name);
    outputFilename = fullfile(outputFolder, eegFilenames(fileCounter).name);
    if(~exist(outputFilename, 'file'))
        eegData = load(filename);
        if(isfield('phaseSlopeTopographies', eegData))
            mriData = load(mriFilename);
            markerTimes = mriData.tvolumemarkers ./ eegFramesPerSecond;
            
            dataSize = size(eegData.phaseSlopeTopographies.estimatedTimeLag);
            
            convolved.timeLag = NaN(length(markerTimes), dataSize(2), dataSize(3));
            convolved.chanlocs = eegData.phaseSlopeTopographies.chanlocs;
            convolved.filename = eegData.phaseSlopeTopographies.filename;
            convolved.details = eegData.phaseSlopeTopographies.details;
            convolved.details.dim1 = 'time';
            convolved.details.dim2 = 'channel';
            convolved.details.dim3 = 'frequencyBand';
            convolved.frequencyLimits = eegData.phaseSlopeTopographies.frequencyLimits;
            
            for i = 1:dataSize(2)
                fprintf('%s: file %d of %d, signal %d of %d\n', char(datetime), fileCounter, length(eegFilenames), i, size(eegData.phaseSlopeTopographies.estimatedTimeLag, 2));
                for j = 1:dataSize(3)
                    signal = eegData.phaseSlopeTopographies.estimatedTimeLag(:,i,j);
                    sampleRate = 10;
                    %can change sample rate to be read from file it
                    %after sampleRates are recorded correctly
                    %sampleRate = data.phaseSlopeTopographies.details.sampleRate;
                    [convolvedSignal, convolvedTimes] = convolveHrf(signal, sampleRate);
                    markerSignals = interp1(convolvedTimes, convolvedSignal, markerTimes);
                    convolved.timeLag(:, i, j) = markerSignals;
                end
            end
            save(outputFilename, 'convolved', '-v7.3');
        end
    end
    
    fileCounter = fileCounter + 1;
end

