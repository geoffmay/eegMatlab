function locations = fileLocations()

locations.rootTemp = 'C:\Users\Neuro\Documents\temp';
if(~exist(locations.rootTemp, 'file'))
    error('temp folder does not exist: %s', locations.rootTemp);
end
locations.eegFourier = fullfile(locations.rootTemp, 'eegFourier');
locations.convolvedEeg = fullfile(locations.rootTemp, 'convolvedEeg');
locations.consolidatedEegFourier = fullfile(locations.rootTemp, 'consolidatedEegFourier');
locations.brainVision = 'C:\Vision\Raw Files\Geoff EEG test\history';
locations.ghermanFolder = 'C:\Users\Neuro\Downloads\GhermanPhiliastides';

fields = fieldnames(locations);
for i = 1:length(fields)
    folder = locations.(fields{i});
    if(~exist(folder))
        mkdir(folder);
    end
end

ghermanFolders = dir(locations.ghermanFolder);
subInd = find(cellfun(@length, strfind({ghermanFolders.name}, 'sub-')) > 0);
counter = 1;
for i = 1:length(subInd)
    nickname = ghermanFolders(subInd(i)).name;
    subFolder = fullfile(locations.ghermanFolder, nickname);
    eegFolder = fullfile(subFolder, 'EEG');
    files = dir(eegFolder);
    eegInd = find(cellfun(@length, strfind({files.name}, 'EEG_data')) > 0);
    for j = 1:length(eegInd)
        runNickname = files(eegInd(j)).name;
        eegPath = fullfile(eegFolder, runNickname);
        if(~strcmp(eegPath(end-4:end), '.html'))
            ghermanFiles{counter} = eegPath;
            counter = counter + 1;
        end
    end
end
locations.ghermanFiles = ghermanFiles;
