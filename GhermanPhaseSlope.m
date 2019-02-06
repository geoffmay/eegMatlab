rootFolder = 'C:\Users\Neuro\Documents\MATLAB\data\GhermanPhilastides\1.0.1\';
outputFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\phaseSlope';


subFolders = dir(rootFolder);
subFolders(cellfun(@length, strfind({subFolders.name}, 'sub-')) == 0) = [];

for (folderCounter = 1:length(subFolders))
    subFolder = fullfile(rootFolder, subFolders(folderCounter).name, 'EEG');
    files = dir(subFolder);
    files(cellfun(@length, strfind({files.name}, 'EEG_data')) == 0) = [];
    for fileCounter = 1:length(files)
        filename = fullfile(subFolder, files(fileCounter).name);
        if(length(strfind(filename, '.html')) == 0)
            outputFilename = fullfile(outputFolder, files(fileCounter).name);
            if(~exist(outputFilename, 'file'))
                placeholder = sprintf('started on %s', char(datetime));
                save(outputFilename, 'placeholder');
                eeg = loadGhermanEeg(filename);
                tic;
                phaseSlopeTopographies = phaseSlopeTopography(eeg);
                elapsed = toc;
                save(outputFilename, 'phaseSlopeTopographies', 'elapsed', '-v7.3');
            end
        end
    end
end

