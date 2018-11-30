dataFolder = '/home/data/EEG/data/boldEeg';
processedFolder = '/home/data/EEG/processed/boldEeg/convolved';



files = dir(dataFolder);
keep = cellfun(@length, strfind({files.name}, '-edf.edf')) > 0;
files = files(keep);

for i = 1:length(files)
    
    inputPath = fullfile(dataFolder, files(i).name);
    outputPath = fullfile(processedFolder, [files(i).name, '.mat']);
    fourierPath = s1brainVisionFourierMatrix(inputPath);
    
    %todo: convolve results (now separate files) with hrf
    %   fourierEeg = load(fourierPath);
    if(~exist(outputPath, 'file'))
        placeHolder = ['started on ', char(datetime)];
        save(outputPath, 'placeHolder');
        %load MRI volume markers
        markerPath = strrep(inputPath, '-edf.edf', '_Pulse Artifact Correction.Markers');
        markerFrames = s6loadVolumeMarkers(markerPath);
        %drop volume markers whose hrf could be contaminated by uncorrected
        %gradient artifact
        [eeg, dropFrames] = s2loadBrainvisionEdf(inputPath);
        dropIndices = find(dropFrames);
        keepMarker = setdiff(markerFrames, dropIndices);
        markerTimes = keepMarker ./ eeg.srate;
        %convolve
        convolvedFourier = s5convolveAtTimepoints(fourierPath, markerTimes);
        convolvedFourier.filename = inputPath;
        convolvedFourier.timeComputed = datetime;
        
        save(outputPath, 'convolvedFourier', '-v7.3');
    end
end