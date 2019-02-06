function marker = convolveEegAtVolumeMarkers(edfFilename, eeg)

markerFolder = 'C:\Vision\Raw Files\Geoff EEG test\history\';
% inputFolder = 'C:\Vision\Raw Files\Geoff EEG test\export\fourier';
% outputFolder = 'C:\Vision\Raw Files\Geoff EEG test\export\fourierHrf';
% inputSubFolder = fullfile(inputFolder, edfFilename);
% convolvedSubFolder = fullfile(outputFolder, edfFilename);
% mkdir(convolvedSubFolder);
% inputFiles = dir(inputSubFolder);

markerFilename = strrep(edfFilename, '-edf.edf', '_Pulse Artifact Correction.Markers');
markers = loadVolumeMarkers(fullfile(markerFolder, markerFilename));
downsampleRate = 5;
dMarkers = markers / downsampleRate;

% for fileCounter = 3:length(inputFiles)
%     inputPath = fullfile(inputSubFolder, inputFiles(fileCounter).name);
%     outputPath = fullfile(convolvedSubFolder, inputFiles(fileCounter).name);
%     if(~exist(outputPath, 'file'))
%         placeholder = sprintf('started on %s', char(datetime));
% save(outputPath, 'placeholder');
% fprintf('convolving fourier measure %d of %d (%s)\n', fileCounter, length(inputFiles), char(datetime));
% clear eegFourier;
% data = load(inputPath);
marker.signal = NaN(size(eeg.data,1), length(markers));
for channelNumber = 1:size(eeg.data,1)
    fprintf('convolving channel %d of %d (%s)\n', channelNumber, size(eeg.data,1), char(datetime));
    sig = eeg.data(channelNumber, :);
    
    dSignal = downsample(sig, downsampleRate);
    dRate = eeg.srate / downsampleRate;
    
    dHrf = convolveHrf(dSignal, dRate);
    for i = 1:length(dMarkers)
        targetFrame = markers(i);
        targetTime = targetFrame / eeg.srate;
        pieceFrame = min(find(eeg.times >= targetTime));
        
        marker.timeSeconds = eeg.times(pieceFrame);
        dFrame = pieceFrame / downsampleRate;
        marker.signal(channelNumber, i) = interp1(1:length(dHrf), dHrf, dFrame);
%         if(i == 1)
%             volumeValues = repmat(marker, [1, length(markers)]);
%         end
%         volumeValues(channelNumber, i) = marker;
    end
%     if(length(strfind(piece.label, '-')) > 0)
%         eegFourier.label = [piece.label, ' coherence'];
%     else
%         eegFourier.label = [piece.label, ' power'];
%     end
%     eegFourier.freqInfo = piece.freqInfo;
%     eegFourier.convolvedValues = volumeValues;
end

%         save(outputPath, 'eegFourier', '-v7.3');
%     end
% end
% 
% convolvedFiles = dir(convolvedSubFolder);
% 
% for i = 3:length(convolvedFiles)
%     fprintf('consolidating finished file %d of %d\n', i, length(convolvedFiles));
%     ind = i-2;
%     data = load(fullfile(convolvedSubFolder, convolvedFiles(i).name));
%     label = sprintf('%s %dHz-%dHz', data.eegFourier.label, data.eegFourier.freqInfo.lowFreq, data.eegFourier.freqInfo.highFreq);
%     if(i == 1)
%         consolidated.freqInfo.samplesPerSecond = data.eegFourier.freqInfo.samplesPerSecond;
%         consolidated.freqInfo.coherenceMemoryDurationSeconds = data.eegFourier.freqInfo.coherenceMemoryDurationSeconds;
%         consolidated.freqInfo.fftWindowDurationSeconds = data.eegFourier.freqInfo.fftWindowDurationSeconds;
%         consolidated.labels = cell(1, length(convolvedFiles)-2);
%         consolidated.timeSeconds = [data.eegFourier.convolvedValues.timeSeconds]';
%         consolidated.signal = NaN(length(consolidated.timeSeconds), length(convolvedFiles)-2);
%     end
%     consolidated.labels{ind} = data.eegFourier.label;
%     sig = [data.eegFourier.convolvedValues.signal]';
%     consolidated.signal(:, ind) = sig;
% end






