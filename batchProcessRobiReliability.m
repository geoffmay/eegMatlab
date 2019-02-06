
reliabilityFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\Robi\reliability';
fourierFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\Robi\fourier';
files = getRobiDataFiles;
isBaseline = cellfun(@length, strfind(files, 'baseline')) > 0;
isExit = cellfun(@length, strfind(files, 'outcome')) > 0;
files = files(isBaseline | isExit);
isOpen = cellfun(@length, strfind(files, 'open')) > 0;
isClosed = cellfun(@length, strfind(files, 'closed')) > 0;
files = [files(isOpen), files(isClosed)];

for i = 1:length(files)
    slashes = strfind(files{i}, '\');
    file = files{i};
    outFile = file(slashes(end-2)+1:slashes(end)-1);
    outFile = strrep(outFile, '\', '_');
    outputPath = fullfile(reliabilityFolder, [outFile, '.mat']);
    fourierPath = fullfile(fourierFolder, [outFile '.mat']);
    clear surfCoh;
    if(~exist(fourierPath, 'file'))
        eeg = loadRobiEeg(files{i});
        surfCoh = deriveRobiCoherenceMatrix(eeg);
        save(fourierPath, 'surfCoh', '-v7.3');
    end
    if(~exist(outputPath, 'file'))
        clear surfCoh;
        load(fourierPath);
        placeHolder = sprintf('started on %s', char(datetime));
        save(outputPath, 'placeHolder');
        frameStep = 128 * 5;
        j = 1;
        frameCount = frameStep * j;
        maxFrames = length(surfCoh.timePoints) * .4;
        while(frameCount < maxFrames)
            [ summary ] = asymptoteCoherenceReliability2( surfCoh, frameCount );
            summaries(j) = summary;
            j = j + 1;
            frameCount = frameStep * j;
        end
        if(exist('summaries', 'var'))
            save(outputPath, 'summaries', '-v7.3');
        end
        clear summaries;
    end
end