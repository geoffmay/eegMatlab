function folder = saveSplitCohFile(edfFilename, cohInfo)

locations = fileLocations
folder = fullfile(locations.eegFourier, edfFilename);
if(~exist(folder, 'file'))
    mkdir(folder);
    for i = 1:length(cohInfo.coh)
        fprintf('saving coherence measure %d of %d\n', i, length(cohInfo.coh));
        label = cohInfo.coh(i).label;
        thisCoh = cohInfo.coh(i).coherence;
        clear piece;
        for j = 1:size(thisCoh,2)
            piece.signal = thisCoh(:,j);
            outputFilename = sprintf('coherence_%s_%d.mat', label, j);
            piece.freqInfo.lowFrequency = cohInfo.freqInfo.lowFrequencies(j);
            piece.freqInfo.highFrequency = cohInfo.freqInfo.highFrequencies(j);
            piece.label = sprintf('coherence_%s_%d-%d_Hz', label, piece.freqInfo.lowFrequency, piece.freqInfo.highFrequency);
            piece.freqInfo.coherenceSampleRateHz = cohInfo.freqInfo.coherenceSampleRateHz;
            piece.freqInfo.coherenceMemoryDurationSeconds = cohInfo.freqInfo.coherenceMemoryDurationSeconds;
            piece.freqInfo.fftWindowDurationSeconds = cohInfo.freqInfo.fftWindowDurationSeconds;
            piece.times = cohInfo.x';
            outputPath = fullfile(folder, outputFilename);
            save(outputPath, 'piece', '-v7.3');
        end
    end
    for i = 1:length(cohInfo.channels)
        fprintf('saving power measure %d of %d\n', i, length(cohInfo.channels));
        label = cohInfo.channels(i).label;
        thisCoh = cohInfo.channels(i).absolutePower;
        for j = 1:size(thisCoh,2)
            piece.signal = thisCoh(:,j);
            outputFilename = sprintf('absolutePower_%s_%d', label, j);
            piece.freqInfo.lowFrequency = cohInfo.freqInfo.lowFrequencies(j);
            piece.freqInfo.highFrequency = cohInfo.freqInfo.highFrequencies(j);
            piece.label = sprintf('absolutePower_%s_%d-%d_Hz', label, piece.freqInfo.lowFrequency, piece.freqInfo.highFrequency);
            piece.freqInfo.coherenceSampleRateHz = cohInfo.freqInfo.coherenceSampleRateHz;
            piece.freqInfo.coherenceMemoryDurationSeconds = cohInfo.freqInfo.coherenceMemoryDurationSeconds;
            piece.freqInfo.fftWindowDurationSeconds = cohInfo.freqInfo.fftWindowDurationSeconds;
            piece.times = cohInfo.x';
            outputPath = fullfile(folder, outputFilename);
            save(outputPath, 'piece', '-v7.3');
        end
    end
end