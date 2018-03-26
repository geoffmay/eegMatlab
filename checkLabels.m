function [ output_args ] = checkLabels( input_args )
%CHECKLABELS Summary of this function goes here
%   Detailed explanation goes here

files = getAllBdfs;
allFiles = [files.Auditory1, files.Auditory2, files.Flanker1, files.Flanker2, files.PTSD];

longCheck = false;
if(longCheck)
    refChan = int32(32);
    summaries = [];
    for i = 1:length(allFiles)
        path = allFiles{i};
        EEG = pop_readbdf(path, {}, 43, refChan, false);
        summary.filename = allFiles{i};
        summary.labels = [EEG.chanlocs.labels];
        summary.celledLabels = {EEG.chanlocs.labels};
        if(length(summaries) == 0)
            summaries = summary;
        else
            summaries(end+1) = summary;
        end
    end
    
    table = tabulate({summaries.labels});
    goodIndices = find(strcmp(table{1,1}, {summaries.labels}));
    badIndices = find(strcmp(table{2,1}, {summaries.labels}));
    goodFilenames = {summaries(goodIndices).filename};
    badFilenames = {summaries(badIndices).filename};
else
    EEG.nbchan = 43;
end

badFilenames = {...
    '/home/data/EEG/data/Auditory1/VM112.1.Tones.bdf',...
    '/home/data/EEG/data/Auditory1/VM113.1.Tones.bdf',...
    '/home/data/EEG/data/Auditory1/VM114.1.Tones.bdf',...
    '/home/data/EEG/data/Auditory2/VM104.2.Tones.bdf',...
    '/home/data/EEG/data/Auditory2/VM108.2.Tones.bdf',...
    '/home/data/EEG/data/Auditory2/VM110.2.Tones.bdf',...
    '/home/data/EEG/data/Auditory2/VM111.2.Tones.bdf',...
    '/home/data/EEG/data/Auditory2/VM112.2.Tones.bdf',...
    '/home/data/EEG/data/Flanker1/VM112.1.ANT.bdf',...
    '/home/data/EEG/data/Flanker1/VM113.1.ANT.bdf',...
    '/home/data/EEG/data/Flanker1/VM114.1.ANT.bdf',...
    '/home/data/EEG/data/Flanker2/VM104.2.ANT.bdf',...
    '/home/data/EEG/data/Flanker2/VM108.2.ANT.bdf',...
    '/home/data/EEG/data/Flanker2/VM110.2.ANT.bdf',...
    '/home/data/EEG/data/Flanker2/VM111.2.ANT.bdf',...
    '/home/data/EEG/data/Flanker2/VM112.2.ANT.bdf'};


%path = allFiles{goodIndices(1)};
path = allFiles{1};
fileId = fopen(path);
goodFile = fread(fileId);
fclose(fileId);
string = char(goodFile)';
labelStart = strfind(string, 'Fp1');
labelEnd = labelStart + EEG.nbchan * 256 - 1;
channelInfo = string(labelStart:labelEnd);

% EEG = pop_readbdf(path, {}, 43, refChan, false);
% goodChannels = {EEG.chanlocs.labels};

for i = 1:length(badFilenames)
    path = badFilenames{i};
    fileId = fopen(path);
    file = fread(fileId);
    fclose(fileId);
    for j = labelStart:labelEnd
        file(j) = goodFile(j);
    end
%     string = char(file)';
%     badLabel = string(labelStart:labelEnd);
    movefile(path, strcat(path, '(old)'));
    
    fileId = fopen(path, 'w');
    fwrite(fileId, file);
    fclose(fileId);

    
%     EEG = pop_readbdf(path, {}, 43, refChan, false);
%     EEG.chanlocs.labels = goodChannels;
%     pop_writebdf();
end

end

