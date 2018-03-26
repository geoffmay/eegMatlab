files = getRobiDataFiles;
base = cellfun(@length, strfind(files, 'baseline eyes open')) > 0;
exit = cellfun(@length, strfind(files, 'outcome eyes open')) > 0;

outputFolder = '/home/data/EEG/processed/Robi/edf';

toProcess = find(base | exit);

EEG.nbchan = 34;
[~, EEG.chanlocs] = antChannelLocs;
EEG.srate = 2048;
for fileCounter = 1:length(toProcess)
    fprintf('\n%s: file %d of %d', char(datetime), fileCounter, length(toProcess));
    filename = files{toProcess(fileCounter)};
    data = loadRobiDataFile(filename);
    EEG.data = data';
    EEG.pnts = size(EEG.data, 2);
    EEG.times = (1:EEG.pnts) ./ EEG.srate;
    EEG.ref = 'linked mastoids';
    EEG.event = [];
    
    %output file
    slashes = strfind(filename, '/');
    name = filename(slashes(end-2)+1:slashes(end)-1);
    name = strrep(name, '/', ' ');
    outputPath = fullfile(outputFolder, sprintf('%s.edf', name));
    pop_writeeeg(EEG, outputPath, 'TYPE', 'EDF');
end
