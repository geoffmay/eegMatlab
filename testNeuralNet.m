eegInputSeconds = 6;

inputFolder = 'C:\Vision\Raw Files\Geoff EEG test\history';
inputFiles = dir(inputFolder);
edfFiles = inputFiles(cellfun(@length, strfind({inputFiles.name}, '.edf'))>0);
edfFiles(cellfun(@length, strfind({edfFiles.name}, '.edfcoh.mat'))>0) = [];
edfFiles(cellfun(@length, strfind({edfFiles.name}, '.edfphase.mat'))>0) = [];

movieNetworkTimecourse = textread('C:\Vision\Raw Files\Geoff EEG test\mri\sub-GeoffEEGTest_Movie_network_timecourses.txt');
restNetworkTimecourse = textread('C:\Vision\Raw Files\Geoff EEG test\mri\sub-GeoffEEGTest_RSFC_network_timecourses.txt');
movieNetworkIndices = movieNetworkTimecourse(1,:);
restNetworkIndices = restNetworkTimecourse(1,:);
movieNetworkTimecourse = movieNetworkTimecourse(2:end,:);
restNetworkTimecourse = restNetworkTimecourse(2:end,:);
restNetworkCounter = 1;
movieNetworkCounter = 1;
metadata = readCsv('C:\Vision\Raw Files\Geoff EEG test\mri\sub-GeoffEEGTest_metadata.csv');
metadataRow = 1;
totalDataCounter = 1;

for fileCounter = 1:length(edfFiles)
    fprintf('\nfile %d of %d', fileCounter, length(edfFiles));
    volumeCounter = 1;
    runVolumeCounter = 1;
    inputPath = fullfile(inputFolder, edfFiles(fileCounter).name);
    markerPath = strrep(inputPath, '-edf.edf', '_Pulse Artifact Correction.Markers');
    if(~exist(markerPath, 'file'))
        error('unable to find marker file %s', markerPath);
    end
    eeg = loadBrainvisionEdf(inputPath);
    %todo: sanity check eeg volume event times against MRI run times.
    volumeMarkerFrames{fileCounter} = loadVolumeMarkers(markerPath);
    volumeMarkerTimes = volumeMarkerFrames{fileCounter} / eeg.srate;
    while(volumeCounter <= length(volumeMarkerTimes))
%         fprintf('\nfile %d of %d, volume %d of %d', fileCounter, length(edfFiles), volumeCounter, length(volumeMarkerTimes));
        %choose the mri network activation data based on the task
        task = metadata{metadataRow, 'task'};
        task = task{1};
        if(strcmp(task, 'RSFC'))
            networkActivations(totalDataCounter, :) = restNetworkTimecourse(restNetworkCounter, :);
            restNetworkCounter = restNetworkCounter + 1;
        elseif(strcmp(task, 'Movie'))
            networkActivations(totalDataCounter, :) = movieNetworkTimecourse(movieNetworkCounter, :);
            movieNetworkCounter = movieNetworkCounter + 1;
        else
            error('unhandled task %s', task);
        end
        
        %get EEG data
        eegEndFrame = volumeMarkerFrames{fileCounter}(volumeCounter);
        eegStartFrame = eegEndFrame - eegInputSeconds * eeg.srate + 1;
        eegs{totalDataCounter} = eeg.data(:, eegStartFrame:eegEndFrame);
        
        %increment counters
        runVolumeCounter = runVolumeCounter + 1;
        if(runVolumeCounter > metadata{metadataRow, 'volumeCount'})
            runVolumeCounter = 1;
            metadataRow = metadataRow + 1;
            %for now, skip short runs
            while(metadataRow <= size(metadata,1) && metadata{metadataRow, 'volumeCount'} < 132)
                volumeCounter = volumeCounter + metadata{metadataRow, 'volumeCount'};
                metadataRow = metadataRow + 1;
            end
        end
        volumeCounter = volumeCounter + 1;
        totalDataCounter = totalDataCounter + 1;
    end
end

mriOutput = networkActivations';
mriOutput(all(isnan(mriOutput), 2),:) = [];
remove = any(isnan(mriOutput));
mriOutput(:,remove) = [];
eegs(remove) = [];

isTraining = rand(1, size(mriOutput, 2)) >= 0.5;
eegInput = NaN(numel(eegs{1}), length(eegs));
for i = 1:length(eegs)
    eegInput(:,i) = reshape(eegs{i}, [numel(eegs{i}), 1]);
end

trainingEeg = eegInput(:, isTraining);
trainingMri = mriOutput(:, isTraining);
testEeg = eegInput(:, ~isTraining);
testMri = mriOutput(:, ~isTraining);



