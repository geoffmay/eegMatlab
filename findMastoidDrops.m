folder = '/media/eegDrive';
files = dir(folder);
files([files.isdir]) = [];
for(fileCounter = 1:length(files))
    outputFilename = fullfile('/home/data/EEG/processed/Robi/overallCoherence', files(fileCounter).name);
    fprintf('\n%s: %d of %d', char(datetime), fileCounter, length(files));
    if(~exist(outputFilename,'file'))
        filename = fullfile(folder, files(fileCounter).name);
        file = load(filename);
        if(isfield(file, 'channelPairs'))
        for i = 1:length(file.channelPairs)
            if(i == 1)
                allCoherence = mean(file.channelPairs(i).coherence,2);
            else
                allCoherence = allCoherence + mean(file.channelPairs(i).coherence,2);
            end
        end
        allCoherence = allCoherence ./ length(file.channelPairs);
        save(outputFilename, 'allCoherence');
        end
    end
end


