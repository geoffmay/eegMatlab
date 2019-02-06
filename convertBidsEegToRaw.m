rootFolder = 'C:\Users\Neuro\Downloads\GhermanPhiliastides';
subFolders = dir(rootFolder);
subFolders(1:2) = [];
subFolders(~[subFolders.isdir]) = [];
for i = 1:length(subFolders)
    eegFolder = fullfile(rootFolder, subFolders(i).name, 'EEG');
    files = dir(eegFolder);
    files(1:2) = [];
    files = files(cellfun(@length, strfind({files.name}, 'EEG_data')) > 0);
    for j = 1:length(files)
        eegPath = fullfile(rootFolder, subFolders(i).name, 'EEG', files(j).name)
        eventPath = strrep(eegPath, 'EEG_data', 'EEG_events');
        eeg = load(eegPath);
        events = load(eventPath);
        eegOutputPath = strrep(eegPath, 'EEG_data_','');
        eegOutputPath = strrep(eegOutputPath, '.mat','.raw');
        saveRawMatrix(eeg.EEGdata.Y, eegOutputPath);
    end
end