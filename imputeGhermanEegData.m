rootFolder = 'C:\Users\Neuro\Documents\MATLAB\data\GhermanPhilastides\1.0.1\';
outputFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\imputed';

subFolders = dir(rootFolder);
subFolders(cellfun(@length, strfind({subFolders.name}, 'sub-')) == 0) = [];
if(~exist(outputFolder, 'file'))
    mkdir(outputFolder);
end

chanlocs = readlocs('C:\Users\Neuro\Documents\MATLAB\data\GhermanPhilastides\1.0.1\additional_files_electrode_info.elp');
for i = 1:length(chanlocs)
    chanlabels{i} = chanlocs(i).labels;
end

%load the data
fileCounter = 1;
sampleCounter = 0;
for (folderCounter = 1:length(subFolders))
    subFolder = fullfile(rootFolder, subFolders(folderCounter).name, 'EEG');
    files = dir(subFolder);
    files(cellfun(@length, strfind({files.name}, 'EEG_data')) == 0) = [];
    for fileCounter = 1:length(files)
        fprintf('loading folder %d file %d', folderCounter, fileCounter);
        filename = fullfile(subFolder, files(fileCounter).name);
        %         outputFilename = fullfile(outputFolder, files(fileCounter).name);
        %         if(~exist(outputFilename, 'file'))
        %placeholder = sprintf('started on %s', char(datetime));
        %         save(outputFilename, 'placeholder');
        if(length(strfind(filename, '.html')) == 0)
            eeg = loadGhermanEeg(filename);
            data = eeg.data';
            data = downsample(data, 100);
            firstFrame = sampleCounter+1;
            lastFrame = sampleCounter + size(data,1);
            bigEeg(firstFrame:lastFrame, :) = NaN;
            for i = 1:eeg.nbchan
                ind = find(strcmp(chanlabels, eeg.chanlocs(i).labels));
                bigEeg(firstFrame:lastFrame, ind) = data(:, i); %eeg.data(i, :)';
            end
            sampleCounter = lastFrame;
            lastFrames(fileCounter) = lastFrame;
            fileCounter = fileCounter + 1;
        end
        %             tic;
        %             phaseSlopeTopographies = phaseSlopeTopography(eeg);
        %             elapsed = toc;
        %             save(outputFilename, 'phaseSlopeTopographies', 'elapsed', '-v7.3');
    end
end

save('C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\eegSrate100.mat', 'bigEeg', '-v7.3');

[output, pcaResult] = imputeUsingPpca(bigEeg);
save('C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\eegSrate100.mat', 'output', 'pcaResult', '-v7.3');
