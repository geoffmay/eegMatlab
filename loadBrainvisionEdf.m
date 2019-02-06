function [eeg, dropFrames] = loadBrainvisionEdf(filename, secondsSweep)
%LOADBRAINVISIONEDF Summary of this function goes here
%   Detailed explanation goes here
eeg = eeg_emptyset;
eeg.data = ft_read_data(filename);
header = ft_read_header(filename);
eeg.srate = header.Fs;
eeg.nbchan = header.nChans;
eeg.pnts = size(eeg.data, 2);
for i = 1:length(header.label)
    eeg.chanlocs(i).labels = header.label{i};
end
eeg = setStandardLocations(eeg);
eeg.times = (1:eeg.pnts) ./ eeg.srate;
[~, file, ext] = fileparts(filename);
eeg.filename = [file ext];
eeg.filepath = filename;

diffSum = zeros(1, size(eeg.data, 2) - 1);
for i = 1:size(eeg.data,1)
    d = diff(eeg.data(1,:));
    diffSum = diffSum + abs(d);
end
diffSum = diffSum ./ size(eeg.data,1);

threshold = 100;
supraThreshold = diffSum > threshold;
if(~exist('secondsSweep', 'var'))
    secondsSweep = [1.5, 7.5];
end
frameSweep = secondsSweep .* eeg.srate;
dropFrames = supraThreshold;

for i = 1:length(supraThreshold)
    if(supraThreshold(i) > 0)
        dropFrames((i-frameSweep(1)):(i+frameSweep(2))) = 1;
    end
end


end

