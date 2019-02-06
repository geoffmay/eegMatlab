

filename = 'C:\Users\Neuro\Downloads\TMS eeg\brainvision\patient001-edf.raw';
sampleRate = 5000;
stimRate = 18;

rereference = false;
channelCount = 64;
file = dir(filename);
fileLength = file.bytes / 8;
sampleCount = fileLength / channelCount;
if(floor(sampleCount) ~= sampleCount)
    sampleCount = floor(sampleCount);
    fileLength = sampleCount * channelCount;
end
[folder, file, ext] = fileparts(filename);
fileTimestamp = str2num(file);
multiplier = 1;
fileId = fopen(filename);
contents = fread(fileId, fileLength, 'double');
fclose(fileId);

%truncate if necessary
sampleCount = fileLength / channelCount;
if(sampleCount ~= floor(sampleCount))
    sampleCount = floor(sampleCount);
    fileLength = sampleCount * channelCount;
    contents = contents(1:fileLength);
end
labels = antChannelLocs;

data = reshape(contents, fileLength / channelCount, channelCount);
%clear contents;

if(rereference)
    cpzIndex = find(strcmp(labels,'CPz'));
    m1Index = find(strcmp(labels,'M1'));
    m2Index = find(strcmp(labels,'M2'));
    
    for i=1:size(data,1)
        sample = data(i,:);
        sample = sample .* multiplier;
        linkedMastoids = (sample(m1Index) + sample(m2Index)) / 2;
        %     avg = mean(sample(1:33));
        newSample = sample - linkedMastoids;
        %     avgSample = sample - avg;
        data(i,1:33) = newSample(1:33);
    end
end

channelLabels = brainvisionTmsChannelLabels;

unused1Index = find(strcmp(channelLabels, 'unused1'));
unused2Index = find(strcmp(channelLabels, 'unused2'));
data(:, [unused1Index, unused2Index]) = [];
channelLabels([unused1Index, unused2Index]) = [];
meanData = mean(abs(data),2);
x = (1:size(data,1)) ./ sampleRate;

%plot
plotAll = false;
if(plotAll)
    xMin = 653;
    xMax = 685;
    frames = 1:size(data,1);
else
    xMin = 653;
    xMax = 685;
    frames = (xMin*sampleRate):(xMax*sampleRate);
end
close all;
plotIndex = [1:2 6];
toPlot = data(frames,plotIndex);
plotLabels = labels(plotIndex);
plot(x(frames), toPlot);
legend(plotLabels);
xlabel('time (seconds)');
ylabel('voltage');


%cut tms pieces
threshold = 600;
pulseTimes = meanData > threshold;
trimSize = ceil(sampleRate / stimRate);
paraPulse = false(size(pulseTimes));
for i = 1:length(pulseTimes)
    if(pulseTimes(i))
        paraPulse(i-trimSize:i+trimSize) = true;
    end
end
