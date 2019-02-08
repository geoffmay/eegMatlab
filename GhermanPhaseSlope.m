rootFolder = 'C:\Users\Neuro\Documents\MATLAB\data\GhermanPhilastides\1.0.1\';
phaseSlopeFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\phaseSlope';

%compute
subFolders = dir(rootFolder);
subFolders(cellfun(@length, strfind({subFolders.name}, 'sub-')) == 0) = [];
for (folderCounter = 1:length(subFolders))
    subFolder = fullfile(rootFolder, subFolders(folderCounter).name, 'EEG');
    files = dir(subFolder);
    files(cellfun(@length, strfind({files.name}, 'EEG_data')) == 0) = [];
    for fileCounter = 1:length(files)
        filename = fullfile(subFolder, files(fileCounter).name);
        if(length(strfind(filename, '.html')) == 0)
            outputFilename = fullfile(phaseSlopeFolder, files(fileCounter).name);
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

%consolidate
phaseFiles = dir(phaseSlopeFolder);
phaseFiles([phaseFiles.isdir]) = [];
allChans = readlocs('C:\Users\Neuro\Documents\MATLAB\data\GhermanPhilastides\1.0.1\additional_files_electrode_info.elp');
allLabels = {allChans.labels};
startIndex = 1;
for i = 1:length(phaseFiles)
    fprintf('');
    filePath = fullfile(phaseSlopeFolder, phaseFiles(i).name);
    a = load(filePath);
    b = a.phaseSlopeTopographies.estimatedTimeLag(:, :, 1);
    endIndex = startIndex + size(b, 1) - 1;
    allTopographies.data(startIndex:endIndex, :) = NaN(size(b,1), length(allLabels));
    for j = 1:size(b, 2)
        labels = {a.phaseSlopeTopographies.chanlocs.labels};
        lab = labels{j};
        if(strcmp(lab(1:2), 'Ch'))
            ind = str2num(lab(3:end));
        else
            ind = find(strcmp(allLabels, lab));
        end
        allTopographies.data(startIndex:endIndex, ind) = b(:, j);
    end
    %     allTopographies.data(startIndex:endIndex, :) = b;
    allTopographies.fileNumber(startIndex:endIndex) = i;
    allTopographies.fileNames{i} = filePath;
    startIndex = endIndex + 1;
end

%impute
[allTopographies.data, pcaResult] = imputeUsingPpca(allTopographies.data);



