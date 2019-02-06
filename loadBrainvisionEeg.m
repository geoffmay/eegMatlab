channelLabels = brainvisionTmsChannelLabels;

filename = 'C:\Users\Neuro\Downloads\TMS eeg\brainvision\patient001-edf.raw';
sampleRate = 5000;
stimRate = 18;

%load data
data = loadBrainvisionDataFile(filename);
unused1Index = find(strcmp(channelLabels, 'unused1'));
unused2Index = find(strcmp(channelLabels, 'unused2'));
data(:, [unused1Index, unused2Index]) = [];
channelLabels([unused1Index, unused2Index]) = [];
meanData = mean(abs(data),2);
x = (1:size(data,1)) ./ sampleRate;

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



%perform ica
[ica.weights,ica.sphere,ica.compvars,ica.bias,ica.signs,ica.lrates,ica.activations] = runica(data');

%look at one channel
displayLabel = 'F6';
index = find(strcmp(channelLabels, displayLabel));
seconds = 393; %manually panning here for now
close all
dat = data(:,index);
plot(x, dat);
pan xon;
zoom xon;
title(sprintf('EEG before, during and after TMS: channel %s', displayLabel));
ylabel('voltage uV');
xlabel('time, seconds');



save('C:\Users\Neuro\Documents\MATLAB\processed\workspace.mat', '-v7.4');