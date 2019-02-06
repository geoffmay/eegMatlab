function [EEG] = loadBrainvisionBinary(filename)
%LOADBRAINVISIONBINARY Summary of this function goes here
%   Detailed explanation goes here

%parse header
if(strcmp(lower(filename(end-3:end)),'.dat'))
    headerName = [filename(1:end-3) 'vhdr'];
    if(~exist(headerName, 'file'))
        error(sprintf('header file does not exist: %s', headerName));
    end
    headerFileId = fopen(headerName);
    EEG.header = fscanf(headerFileId, '%c');
    fclose(headerFileId);
    EEG.chanCount = parseLineValue(EEG.header, 'NumberOfChannels=');
    EEG.dataPoints = parseLineValue(EEG.header, 'DataPoints=');
    EEG.sampleIntervalNanoseconds = parseLineValue(EEG.header, 'SamplingInterval=');
    EEG.sampleRate = 1000000 / EEG.sampleIntervalNanoseconds;
    EEG.dataOrientation = parseLineValue(EEG.header, 'DataOrientation=');
    EEG.dataType = parseLineValue(EEG.header, 'DataType=');
    EEG.dataFormat = parseLineValue(EEG.header, 'DataFormat=');
    EEG.binaryFormat = parseLineValue(EEG.header, 'BinaryFormat=');
    EEG.channelLabels = cell(1,EEG.chanCount);
    for i = 1:EEG.chanCount
        chan = sprintf('Ch%d=', i);
        chanLine = parseLineValue(EEG.header, chan);
        chanItems = strsplit(chanLine, ',');
        EEG.channelLabels{i} = chanItems{1};
    end
end
EEG.chanlocs = getChanlocs(EEG.channelLabels);

%read data
file = dir(filename);
fileLength = file.bytes / 4;
fileId = fopen(filename);
contents = fread(fileId, fileLength, 'single');
fclose(fileId);
sampleCount = length(contents) / EEG.chanCount;
if(strcmp(EEG.dataOrientation, 'VECTORIZED'))
    EEG.data = reshape(contents, [sampleCount, EEG.chanCount]);
else
    error('unsupported data orientation, need to write code for multiplexed here.');
end
EEG.timepoints = (1:size(EEG.data,1)) ./ EEG.sampleRate;
end

