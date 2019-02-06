function [ channelPairs, x, channels, freqInfo ] = allChannelCoherence3( EEG )
%COHERENCE Computes a timeseries of coherence for a given frequency.
%   Detailed explanation goes here

highRes = 0;
doRelative = 0;
doRatios = 1;

processFrequenciesInGroups = 1;
if(highRes)
  fftWindowDuration = 4;
else
  fftWindowDuration = 1;
end

if(processFrequenciesInGroups)
  
  if(highRes)
    freqDiv = [.25 .5 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22 24 26 28 30 32 36 40 44 48 52 56 60 64];
    lowFreq = freqDiv(1:end-1);
    highFreq = freqDiv(2:end)-(1/fftWindowDuration);
  else
    lowFreq = [1 5 9 13 25];
    highFreq = [4 8 12 24 30];
    %         lowFreq = [1 4 8 12 25];
    %         highFreq = [4 8 12 25 30];
  end
  
  freqCount = length(lowFreq);
  iRange = 1:freqCount;
else
  lowFreq = 1/fftWindowDuration;
  highFreq = 64;
  freqCount = (highFreq - lowFreq + 1) * fftWindowDuration;
  iRange = 2:freqCount+1;
end
sampleRate = EEG.srate;
refreshRate = 128;
coherenceMemoryDuration = 0.5;

freqInfo.lowFrequencies = lowFreq;
freqInfo.highFrequencies = highFreq;
freqInfo.coherenceSampleRateHz = refreshRate;
freqInfo.coherenceMemoryDurationSeconds = coherenceMemoryDuration;
freqInfo.fftWindowDurationSeconds = fftWindowDuration;


debugMode = false;


%initialize some vars
if(~exist('EEG', 'var'))
  error('No data was passed to this function');
end

plotLength=floor(size(EEG.data,2) / (EEG.srate/refreshRate) - refreshRate * fftWindowDuration);
channelPairCount = EEG.nbchan * (EEG.nbchan-1) / 2;
nullPlot = NaN(plotLength,1);

% channelPairs.coherence = NaN(plotLength, freqCount,channelPairCount);
% channels.absolutePower = NaN(plotLength, freqCount,EEG.nbchan);
channel.absolutePower = nullPlot;
for i = 1:length(EEG.chanlocs)
  channel.label = EEG.chanlocs(i).labels;
  channels(i) = channel;
end

%channelPairs.labels = cell(1,length(channels.labels) * (length(channels.labels)-1) / 2);
channelPair.coherence = nullPlot;
counter = 1;
for i = 1:length(channels)
  for j = i+1:length(channels)
    channelPair.label = sprintf('%s-%s', channels(i).label, channels(j).label);
    channelPairs(counter) = channelPair;
    counter = counter + 1;
  end
end

tic;

outputFolder = '/home/data/EEG/processed/Robi/coherence3';
tempFolder = '/home/data/EEG/processed/Robi/temp';
headerFilename = 'ROBI_003_baseline_eyes_open1.coh.txt';
binaryFilename = 'ROBI_003_baseline_eyes_open1.coh';
fileId = fopen(fullfile(outputFolder, headerFilename));
header = fscanf(fileId, '%c');
lines = strfind(header, sprintf('\n'));
output.sRate = sscanf(header(1:(lines(1)-1)), 'sample rate: %f');
output.cohRate = sscanf(header((lines(1)+1):(lines(2)-1)), 'coherence refresh rate: %f');
interp = header((lines(2)+1):(lines(3)-1));
output.interpolation = interp(length('interpolation: ')+1:end);
output.measureCount = sscanf(header((lines(3)+1):(lines(4)-1)), 'measure count: %f');
measures = header((lines(4)+1):(lines(5)-1));
measures = measures(length('measurements (comma separated): ')+1:end);
output.measures = strsplit(measures, ',');
output.data = loadBinaryMatrix(fullfile(outputFolder,binaryFilename), output.measureCount);
wastedColumns = find(all(isnan(output.data),1));
incompleteRows = find(any(isnan(output.data),2));

template = sprintf('sample rate: %%f\ncoherence refresh rate: %%f
measure count: 2325
measurements (comma separated): abs Fp1 1Hz-3Hz,abs Fp1 4Hz-7Hz,abs Fp1 8Hz-11Hz,abs Fp1 12Hz-24Hz,abs Fp1 25Hz-29Hz,abs Fpz 1Hz-3Hz,abs Fpz 4Hz-7Hz,abs Fpz 8Hz-11Hz,abs Fpz 12Hz-24Hz,abs Fpz 25Hz-29Hz,abs Fp2 1Hz-3Hz,abs Fp2 4Hz-7Hz,abs Fp2 8Hz-11Hz,abs Fp2 12Hz-24Hz,abs Fp2 25Hz-29Hz,abs F7 1Hz-3Hz,abs F7 4Hz-7Hz,abs F7 8Hz-11Hz,abs F7 12Hz-24Hz,abs F7 25Hz-29Hz,abs F3 1Hz-3Hz,abs F3 4Hz-7Hz,abs F3 8Hz-11Hz,abs F3 12Hz-24Hz,abs F3 25Hz-29Hz,abs Fz 1Hz-3Hz,abs Fz 4Hz-7Hz,abs Fz 8Hz-11Hz,abs Fz 12Hz-24Hz,abs Fz 25Hz-29Hz,abs F4 1Hz-3Hz,abs F4 4Hz-7Hz,abs F4 8Hz-11Hz,abs F4 12Hz-24Hz,abs F4 25Hz-29Hz,abs F8 1Hz-3Hz,abs F8 4Hz-7Hz,abs F8 8Hz-11Hz,abs F8 12Hz-24Hz,abs F8 25Hz-29Hz,abs FC5 1Hz-3Hz,abs FC5 4Hz-7Hz,abs FC5 8Hz-11Hz,abs FC5 12Hz-24Hz,abs FC5 25Hz-29Hz,abs FC1 1Hz-3Hz,abs FC1 4Hz-7Hz,abs FC1 8Hz-11Hz,abs FC1 12Hz-24Hz,abs FC1 25Hz-29Hz,abs FC2 1Hz-3Hz,abs FC2 4Hz-7Hz,abs FC2 8Hz-11Hz,abs FC2 12Hz-24Hz,abs FC2 25Hz-29Hz,abs FC6 1Hz-3Hz,abs FC6 4Hz-7Hz,abs FC6 8Hz-11Hz,abs FC6 12Hz-24Hz,abs FC6 25Hz-29Hz,abs T7 1Hz-3Hz,abs T7 4Hz-7Hz,abs T7 8Hz-11Hz,abs T7 12Hz-24Hz,abs T7 25Hz-29Hz,abs C3 1Hz-3Hz,abs C3 4Hz-7Hz,abs C3 8Hz-11Hz,abs C3 12Hz-24Hz,abs C3 25Hz-29Hz,abs Cz 1Hz-3Hz,abs Cz 4Hz-7Hz,abs Cz 8Hz-11Hz,abs Cz 12Hz-24Hz,abs Cz 25Hz-29Hz,abs C4 1Hz-3Hz,abs C4 4Hz-7Hz,abs C4 8Hz-11Hz,abs C4 12Hz-24Hz,abs C4 25Hz-29Hz,abs T8 1Hz-3Hz,abs T8 4Hz-7Hz,abs T8 8Hz-11Hz,abs T8 12Hz-24Hz,abs T8 25Hz-29Hz,abs CP5 1Hz-3Hz,abs CP5 4Hz-7Hz,abs CP5 8Hz-11Hz,abs CP5 12Hz-24Hz,abs CP5 25Hz-29Hz,abs CP1 1Hz-3Hz,abs CP1 4Hz-7Hz,abs CP1 8Hz-11Hz,abs CP1 12Hz-24Hz,abs CP1 25Hz-29Hz,abs CP2 1Hz-3Hz,abs CP2 4Hz-7Hz,abs CP2 8Hz-11Hz,abs CP2 12Hz-24Hz,abs CP2 25Hz-29Hz,abs CP6 1Hz-3Hz,abs CP6 4Hz-7Hz,abs CP6 8Hz-11Hz,abs CP6 12Hz-24Hz,abs CP6 25Hz-29Hz,abs P7 1Hz-3Hz,abs P7 4Hz-7Hz,abs P7 8Hz-11Hz,abs P7 12Hz-24Hz,abs P7 25Hz-29Hz,abs P3 1Hz-3Hz,abs P3 4Hz-7Hz,abs P3 8Hz-11Hz,abs P3 12Hz-24Hz,abs P3 25Hz-29Hz,abs Pz 1Hz-3Hz,abs Pz 4Hz-7Hz,abs Pz 8Hz-11Hz,abs Pz 12Hz-24Hz,abs Pz 25Hz-29Hz,abs P4 1Hz-3Hz,abs P4 4Hz-7Hz,abs P4 8Hz-11Hz,abs P4 12Hz-24Hz,abs P4 25Hz-29Hz,abs P8 1Hz-3Hz,abs P8 4Hz-7Hz,abs P8 8Hz-11Hz,abs P8 12Hz-24Hz,abs P8 25Hz-29Hz,abs POz 1Hz-3Hz,abs POz 4Hz-7Hz,abs POz 8Hz-11Hz,abs POz 12Hz-24Hz,abs POz 25Hz-29Hz,abs O1 1Hz-3Hz,abs O1 4Hz-7Hz,abs O1 8Hz-11Hz,abs O1 12Hz-24Hz,abs O1 25Hz-29Hz,abs Oz 1Hz-3Hz,abs Oz 4Hz-7Hz,abs Oz 8Hz-11Hz,abs Oz 12Hz-24Hz,abs Oz 25Hz-29Hz,abs O2 1Hz-3Hz,abs O2 4Hz-7Hz,abs O2 8Hz-11Hz,abs O2 12Hz-24Hz,abs O2 25Hz-29Hz,coh Fp1-Fpz 1Hz-3Hz,coh Fp1-Fpz 4Hz-7Hz,coh Fp1-Fpz 8Hz-11Hz,coh Fp1-Fpz 12Hz-24Hz,coh Fp1-Fpz 25Hz-29Hz,coh Fp1-Fp2 1Hz-3Hz,coh Fp1-Fp2 4Hz-7Hz,coh Fp1-Fp2 8Hz-11Hz,coh Fp1-Fp2 12Hz-24Hz,coh Fp1-Fp2 25Hz-29Hz,coh Fp1-F7 1Hz-3Hz,coh Fp1-F7 4Hz-7Hz,coh Fp1-F7 8Hz-11Hz,coh Fp1-F7 12Hz-24Hz,coh Fp1-F7 25Hz-29Hz,coh Fp1-F3 1Hz-3Hz,coh Fp1-F3 4Hz-7Hz,coh Fp1-F3 8Hz-11Hz,coh Fp1-F3 12Hz-24Hz,coh Fp1-F3 25Hz-29Hz,coh Fp1-Fz 1Hz-3Hz,coh Fp1-Fz 4Hz-7Hz,coh Fp1-Fz 8Hz-11Hz,coh Fp1-Fz 12Hz-24Hz,coh Fp1-Fz 25Hz-29Hz,coh Fp1-F4 1Hz-3Hz,coh Fp1-F4 4Hz-7Hz,coh Fp1-F4 8Hz-11Hz,coh Fp1-F4 12Hz-24Hz,coh Fp1-F4 25Hz-29Hz,coh Fp1-F8 1Hz-3Hz,coh Fp1-F8 4Hz-7Hz,coh Fp1-F8 8Hz-11Hz,coh Fp1-F8 12Hz-24Hz,coh Fp1-F8 25Hz-29Hz,coh Fp1-FC5 1Hz-3Hz,coh Fp1-FC5 4Hz-7Hz,coh Fp1-FC5 8Hz-11Hz,coh Fp1-FC5 12Hz-24Hz,coh Fp1-FC5 25Hz-29Hz,coh Fp1-FC1 1Hz-3Hz,coh Fp1-FC1 4Hz-7Hz,coh Fp1-FC1 8Hz-11Hz,coh Fp1-FC1 12Hz-24Hz,coh Fp1-FC1 25Hz-29Hz,coh Fp1-FC2 1Hz-3Hz,coh Fp1-FC2 4Hz-7Hz,coh Fp1-FC2 8Hz-11Hz,coh Fp1-FC2 12Hz-24Hz,coh Fp1-FC2 25Hz-29Hz,coh Fp1-FC6 1Hz-3Hz,coh Fp1-FC6 4Hz-7Hz,coh Fp1-FC6 8Hz-11Hz,coh Fp1-FC6 12Hz-24Hz,coh Fp1-FC6 25Hz-29Hz,coh Fp1-T7 1Hz-3Hz,coh Fp1-T7 4Hz-7Hz,coh Fp1-T7 8Hz-11Hz,coh Fp1-T7 12Hz-24Hz,coh Fp1-T7 25Hz-29Hz,coh Fp1-C3 1Hz-3Hz,coh Fp1-C3 4Hz-7Hz,coh Fp1-C3 8Hz-11Hz,coh Fp1-C3 12Hz-24Hz,coh Fp1-C3 25Hz-29Hz,coh Fp1-Cz 1Hz-3Hz,coh Fp1-Cz 4Hz-7Hz,coh Fp1-Cz 8Hz-11Hz,coh Fp1-Cz 12Hz-24Hz,coh Fp1-Cz 25Hz-29Hz,coh Fp1-C4 1Hz-3Hz,coh Fp1-C4 4Hz-7Hz,coh Fp1-C4 8Hz-11Hz,coh Fp1-C4 12Hz-24Hz,coh Fp1-C4 25Hz-29Hz,coh Fp1-T8 1Hz-3Hz,coh Fp1-T8 4Hz-7Hz,coh Fp1-T8 8Hz-11Hz,coh Fp1-T8 12Hz-24Hz,coh Fp1-T8 25Hz-29Hz,coh Fp1-CP5 1Hz-3Hz,coh Fp1-CP5 4Hz-7Hz,coh Fp1-CP5 8Hz-11Hz,coh Fp1-CP5 12Hz-24Hz,coh Fp1-CP5 25Hz-29Hz,coh Fp1-CP1 1Hz-3Hz,coh Fp1-CP1 4Hz-7Hz,coh Fp1-CP1 8Hz-11Hz,coh Fp1-CP1 12Hz-24Hz,coh Fp1-CP1 25Hz-29Hz,coh Fp1-CP2 1Hz-3Hz,coh Fp1-CP2 4Hz-7Hz,coh Fp1-CP2 8Hz-11Hz,coh Fp1-CP2 12Hz-24Hz,coh Fp1-CP2 25Hz-29Hz,coh Fp1-CP6 1Hz-3Hz,coh Fp1-CP6 4Hz-7Hz,coh Fp1-CP6 8Hz-11Hz,coh Fp1-CP6 12Hz-24Hz,coh Fp1-CP6 25Hz-29Hz,coh Fp1-P7 1Hz-3Hz,coh Fp1-P7 4Hz-7Hz,coh Fp1-P7 8Hz-11Hz,coh Fp1-P7 12Hz-24Hz,coh Fp1-P7 25Hz-29Hz,coh Fp1-P3 1Hz-3Hz,coh Fp1-P3 4Hz-7Hz,coh Fp1-P3 8Hz-11Hz,coh Fp1-P3 12Hz-24Hz,coh Fp1-P3 25Hz-29Hz,coh Fp1-Pz 1Hz-3Hz,coh Fp1-Pz 4Hz-7Hz,coh Fp1-Pz 8Hz-11Hz,coh Fp1-Pz 12Hz-24Hz,coh Fp1-Pz 25Hz-29Hz,coh Fp1-P4 1Hz-3Hz,coh Fp1-P4 4Hz-7Hz,coh Fp1-P4 8Hz-11Hz,coh Fp1-P4 12Hz-24Hz,coh Fp1-P4 25Hz-29Hz,coh Fp1-P8 1Hz-3Hz,coh Fp1-P8 4Hz-7Hz,coh Fp1-P8 8Hz-11Hz,coh Fp1-P8 12Hz-24Hz,coh Fp1-P8 25Hz-29Hz,coh Fp1-POz 1Hz-3Hz,coh Fp1-POz 4Hz-7Hz,coh Fp1-POz 8Hz-11Hz,coh Fp1-POz 12Hz-24Hz,coh Fp1-POz 25Hz-29Hz,coh Fp1-O1 1Hz-3Hz,coh Fp1-O1 4Hz-7Hz,coh Fp1-O1 8Hz-11Hz,coh Fp1-O1 12Hz-24Hz,coh Fp1-O1 25Hz-29Hz,coh Fp1-Oz 1Hz-3Hz,coh Fp1-Oz 4Hz-7Hz,coh Fp1-Oz 8Hz-11Hz,coh Fp1-Oz 12Hz-24Hz,coh Fp1-Oz 25Hz-29Hz,coh Fp1-O2 1Hz-3Hz,coh Fp1-O2 4Hz-7Hz,coh Fp1-O2 8Hz-11Hz,coh Fp1-O2 12Hz-24Hz,coh Fp1-O2 25Hz-29Hz,coh Fpz-Fp2 1Hz-3Hz,coh Fpz-Fp2 4Hz-7Hz,coh Fpz-Fp2 8Hz-11Hz,coh Fpz-Fp2 12Hz-24Hz,coh Fpz-Fp2 25Hz-29Hz,coh Fpz-F7 1Hz-3Hz,coh Fpz-F7 4Hz-7Hz,coh Fpz-F7 8Hz-11Hz,coh Fpz-F7 12Hz-24Hz,coh Fpz-F7 25Hz-29Hz,coh Fpz-F3 1Hz-3Hz,coh Fpz-F3 4Hz-7Hz,coh Fpz-F3 8Hz-11Hz,coh Fpz-F3 12Hz-24Hz,coh Fpz-F3 25Hz-29Hz,coh Fpz-Fz 1Hz-3Hz,coh Fpz-Fz 4Hz-7Hz,coh Fpz-Fz 8Hz-11Hz,coh Fpz-Fz 12Hz-24Hz,coh Fpz-Fz 25Hz-29Hz,coh Fpz-F4 1Hz-3Hz,coh Fpz-F4 4Hz-7Hz,coh Fpz-F4 8Hz-11Hz,coh Fpz-F4 12Hz-24Hz,coh Fpz-F4 25Hz-29Hz,coh Fpz-F8 1Hz-3Hz,coh Fpz-F8 4Hz-7Hz,coh Fpz-F8 8Hz-11Hz,coh Fpz-F8 12Hz-24Hz,coh Fpz-F8 25Hz-29Hz,coh Fpz-FC5 1Hz-3Hz,coh Fpz-FC5 4Hz-7Hz,coh Fpz-FC5 8Hz-11Hz,coh Fpz-FC5 12Hz-24Hz,coh Fpz-FC5 25Hz-29Hz,coh Fpz-FC1 1Hz-3Hz,coh Fpz-FC1 4Hz-7Hz,coh Fpz-FC1 8Hz-11Hz,coh Fpz-FC1 12Hz-24Hz,coh Fpz-FC1 25Hz-29Hz,coh Fpz-FC2 1Hz-3Hz,coh Fpz-FC2 4Hz-7Hz,coh Fpz-FC2 8Hz-11Hz,coh Fpz-FC2 12Hz-24Hz,coh Fpz-FC2 25Hz-29Hz,coh Fpz-FC6 1Hz-3Hz,coh Fpz-FC6 4Hz-7Hz,coh Fpz-FC6 8Hz-11Hz,coh Fpz-FC6 12Hz-24Hz,coh Fpz-FC6 25Hz-29Hz,coh Fpz-T7 1Hz-3Hz,coh Fpz-T7 4Hz-7Hz,coh Fpz-T7 8Hz-11Hz,coh Fpz-T7 12Hz-24Hz,coh Fpz-T7 25Hz-29Hz,coh Fpz-C3 1Hz-3Hz,coh Fpz-C3 4Hz-7Hz,coh Fpz-C3 8Hz-11Hz,coh Fpz-C3 12Hz-24Hz,coh Fpz-C3 25Hz-29Hz,coh Fpz-Cz 1Hz-3Hz,coh Fpz-Cz 4Hz-7Hz,coh Fpz-Cz 8Hz-11Hz,coh Fpz-Cz 12Hz-24Hz,coh Fpz-Cz 25Hz-29Hz,coh Fpz-C4 1Hz-3Hz,coh Fpz-C4 4Hz-7Hz,coh Fpz-C4 8Hz-11Hz,coh Fpz-C4 12Hz-24Hz,coh Fpz-C4 25Hz-29Hz,coh Fpz-T8 1Hz-3Hz,coh Fpz-T8 4Hz-7Hz,coh Fpz-T8 8Hz-11Hz,coh Fpz-T8 12Hz-24Hz,coh Fpz-T8 25Hz-29Hz,coh Fpz-CP5 1Hz-3Hz,coh Fpz-CP5 4Hz-7Hz,coh Fpz-CP5 8Hz-11Hz,coh Fpz-CP5 12Hz-24Hz,coh Fpz-CP5 25Hz-29Hz,coh Fpz-CP1 1Hz-3Hz,coh Fpz-CP1 4Hz-7Hz,coh Fpz-CP1 8Hz-11Hz,coh Fpz-CP1 12Hz-24Hz,coh Fpz-CP1 25Hz-29Hz,coh Fpz-CP2 1Hz-3Hz,coh Fpz-CP2 4Hz-7Hz,coh Fpz-CP2 8Hz-11Hz,coh Fpz-CP2 12Hz-24Hz,coh Fpz-CP2 25Hz-29Hz,coh Fpz-CP6 1Hz-3Hz,coh Fpz-CP6 4Hz-7Hz,coh Fpz-CP6 8Hz-11Hz,coh Fpz-CP6 12Hz-24Hz,coh Fpz-CP6 25Hz-29Hz,coh Fpz-P7 1Hz-3Hz,coh Fpz-P7 4Hz-7Hz,coh Fpz-P7 8Hz-11Hz,coh Fpz-P7 12Hz-24Hz,coh Fpz-P7 25Hz-29Hz,coh Fpz-P3 1Hz-3Hz,coh Fpz-P3 4Hz-7Hz,coh Fpz-P3 8Hz-11Hz,coh Fpz-P3 12Hz-24Hz,coh Fpz-P3 25Hz-29Hz,coh Fpz-Pz 1Hz-3Hz,coh Fpz-Pz 4Hz-7Hz,coh Fpz-Pz 8Hz-11Hz,coh Fpz-Pz 12Hz-24Hz,coh Fpz-Pz 25Hz-29Hz,coh Fpz-P4 1Hz-3Hz,coh Fpz-P4 4Hz-7Hz,coh Fpz-P4 8Hz-11Hz,coh Fpz-P4 12Hz-24Hz,coh Fpz-P4 25Hz-29Hz,coh Fpz-P8 1Hz-3Hz,coh Fpz-P8 4Hz-7Hz,coh Fpz-P8 8Hz-11Hz,coh Fpz-P8 12Hz-24Hz,coh Fpz-P8 25Hz-29Hz,coh Fpz-POz 1Hz-3Hz,coh Fpz-POz 4Hz-7Hz,coh Fpz-POz 8Hz-11Hz,coh Fpz-POz 12Hz-24Hz,coh Fpz-POz 25Hz-29Hz,coh Fpz-O1 1Hz-3Hz,coh Fpz-O1 4Hz-7Hz,coh Fpz-O1 8Hz-11Hz,coh Fpz-O1 12Hz-24Hz,coh Fpz-O1 25Hz-29Hz,coh Fpz-Oz 1Hz-3Hz,coh Fpz-Oz 4Hz-7Hz,coh Fpz-Oz 8Hz-11Hz,coh Fpz-Oz 12Hz-24Hz,coh Fpz-Oz 25Hz-29Hz,coh Fpz-O2 1Hz-3Hz,coh Fpz-O2 4Hz-7Hz,coh Fpz-O2 8Hz-11Hz,coh Fpz-O2 12Hz-24Hz,coh Fpz-O2 25Hz-29Hz,coh Fp2-F7 1Hz-3Hz,coh Fp2-F7 4Hz-7Hz,coh Fp2-F7 8Hz-11Hz,coh Fp2-F7 12Hz-24Hz,coh Fp2-F7 25Hz-29Hz,coh Fp2-F3 1Hz-3Hz,coh Fp2-F3 4Hz-7Hz,coh Fp2-F3 8Hz-11Hz,coh Fp2-F3 12Hz-24Hz,coh Fp2-F3 25Hz-29Hz,coh Fp2-Fz 1Hz-3Hz,coh Fp2-Fz 4Hz-7Hz,coh Fp2-Fz 8Hz-11Hz,coh Fp2-Fz 12Hz-24Hz,coh Fp2-Fz 25Hz-29Hz,coh Fp2-F4 1Hz-3Hz,coh Fp2-F4 4Hz-7Hz,coh Fp2-F4 8Hz-11Hz,coh Fp2-F4 12Hz-24Hz,coh Fp2-F4 25Hz-29Hz,coh Fp2-F8 1Hz-3Hz,coh Fp2-F8 4Hz-7Hz,coh Fp2-F8 8Hz-11Hz,coh Fp2-F8 12Hz-24Hz,coh Fp2-F8 25Hz-29Hz,coh Fp2-FC5 1Hz-3Hz,coh Fp2-FC5 4Hz-7Hz,coh Fp2-FC5 8Hz-11Hz,coh Fp2-FC5 12Hz-24Hz,coh Fp2-FC5 25Hz-29Hz,coh Fp2-FC1 1Hz-3Hz,coh Fp2-FC1 4Hz-7Hz,coh Fp2-FC1 8Hz-11Hz,coh Fp2-FC1 12Hz-24Hz,coh Fp2-FC1 25Hz-29Hz,coh Fp2-FC2 1Hz-3Hz,coh Fp2-FC2 4Hz-7Hz,coh Fp2-FC2 8Hz-11Hz,coh Fp2-FC2 12Hz-24Hz,coh Fp2-FC2 25Hz-29Hz,coh Fp2-FC6 1Hz-3Hz,coh Fp2-FC6 4Hz-7Hz,coh Fp2-FC6 8Hz-11Hz,coh Fp2-FC6 12Hz-24Hz,coh Fp2-FC6 25Hz-29Hz,coh Fp2-T7 1Hz-3Hz,coh Fp2-T7 4Hz-7Hz,coh Fp2-T7 8Hz-11Hz,coh Fp2-T7 12Hz-24Hz,coh Fp2-T7 25Hz-29Hz,coh Fp2-C3 1Hz-3Hz,coh Fp2-C3 4Hz-7Hz,coh Fp2-C3 8Hz-11Hz,coh Fp2-C3 12Hz-24Hz,coh Fp2-C3 25Hz-29Hz,coh Fp2-Cz 1Hz-3Hz,coh Fp2-Cz 4Hz-7Hz,coh Fp2-Cz 8Hz-11Hz,coh Fp2-Cz 12Hz-24Hz,coh Fp2-Cz 25Hz-29Hz,coh Fp2-C4 1Hz-3Hz,coh Fp2-C4 4Hz-7Hz,coh Fp2-C4 8Hz-11Hz,coh Fp2-C4 12Hz-24Hz,coh Fp2-C4 25Hz-29Hz,coh Fp2-T8 1Hz-3Hz,coh Fp2-T8 4Hz-7Hz,coh Fp2-T8 8Hz-11Hz,coh Fp2-T8 12Hz-24Hz,coh Fp2-T8 25Hz-29Hz,coh Fp2-CP5 1Hz-3Hz,coh Fp2-CP5 4Hz-7Hz,coh Fp2-CP5 8Hz-11Hz,coh Fp2-CP5 12Hz-24Hz,coh Fp2-CP5 25Hz-29Hz,coh Fp2-CP1 1Hz-3Hz,coh Fp2-CP1 4Hz-7Hz,coh Fp2-CP1 8Hz-11Hz,coh Fp2-CP1 12Hz-24Hz,coh Fp2-CP1 25Hz-29Hz,coh Fp2-CP2 1Hz-3Hz,coh Fp2-CP2 4Hz-7Hz,coh Fp2-CP2 8Hz-11Hz,coh Fp2-CP2 12Hz-24Hz,coh Fp2-CP2 25Hz-29Hz,coh Fp2-CP6 1Hz-3Hz,coh Fp2-CP6 4Hz-7Hz,coh Fp2-CP6 8Hz-11Hz,coh Fp2-CP6 12Hz-24Hz,coh Fp2-CP6 25Hz-29Hz,coh Fp2-P7 1Hz-3Hz,coh Fp2-P7 4Hz-7Hz,coh Fp2-P7 8Hz-11Hz,coh Fp2-P7 12Hz-24Hz,coh Fp2-P7 25Hz-29Hz,coh Fp2-P3 1Hz-3Hz,coh Fp2-P3 4Hz-7Hz,coh Fp2-P3 8Hz-11Hz,coh Fp2-P3 12Hz-24Hz,coh Fp2-P3 25Hz-29Hz,coh Fp2-Pz 1Hz-3Hz,coh Fp2-Pz 4Hz-7Hz,coh Fp2-Pz 8Hz-11Hz,coh Fp2-Pz 12Hz-24Hz,coh Fp2-Pz 25Hz-29Hz,coh Fp2-P4 1Hz-3Hz,coh Fp2-P4 4Hz-7Hz,coh Fp2-P4 8Hz-11Hz,coh Fp2-P4 12Hz-24Hz,coh Fp2-P4 25Hz-29Hz,coh Fp2-P8 1Hz-3Hz,coh Fp2-P8 4Hz-7Hz,coh Fp2-P8 8Hz-11Hz,coh Fp2-P8 12Hz-24Hz,coh Fp2-P8 25Hz-29Hz,coh Fp2-POz 1Hz-3Hz,coh Fp2-POz 4Hz-7Hz,coh Fp2-POz 8Hz-11Hz,coh Fp2-POz 12Hz-24Hz,coh Fp2-POz 25Hz-29Hz,coh Fp2-O1 1Hz-3Hz,coh Fp2-O1 4Hz-7Hz,coh Fp2-O1 8Hz-11Hz,coh Fp2-O1 12Hz-24Hz,coh Fp2-O1 25Hz-29Hz,coh Fp2-Oz 1Hz-3Hz,coh Fp2-Oz 4Hz-7Hz,coh Fp2-Oz 8Hz-11Hz,coh Fp2-Oz 12Hz-24Hz,coh Fp2-Oz 25Hz-29Hz,coh Fp2-O2 1Hz-3Hz,coh Fp2-O2 4Hz-7Hz,coh Fp2-O2 8Hz-11Hz,coh Fp2-O2 12Hz-24Hz,coh Fp2-O2 25Hz-29Hz,coh F7-F3 1Hz-3Hz,coh F7-F3 4Hz-7Hz,coh F7-F3 8Hz-11Hz,coh F7-F3 12Hz-24Hz,coh F7-F3 25Hz-29Hz,coh F7-Fz 1Hz-3Hz,coh F7-Fz 4Hz-7Hz,coh F7-Fz 8Hz-11Hz,coh F7-Fz 12Hz-24Hz,coh F7-Fz 25Hz-29Hz,coh F7-F4 1Hz-3Hz,coh F7-F4 4Hz-7Hz,coh F7-F4 8Hz-11Hz,coh F7-F4 12Hz-24Hz,coh F7-F4 25Hz-29Hz,coh F7-F8 1Hz-3Hz,coh F7-F8 4Hz-7Hz,coh F7-F8 8Hz-11Hz,coh F7-F8 12Hz-24Hz,coh F7-F8 25Hz-29Hz,coh F7-FC5 1Hz-3Hz,coh F7-FC5 4Hz-7Hz,coh F7-FC5 8Hz-11Hz,coh F7-FC5 12Hz-24Hz,coh F7-FC5 25Hz-29Hz,coh F7-FC1 1Hz-3Hz,coh F7-FC1 4Hz-7Hz,coh F7-FC1 8Hz-11Hz,coh F7-FC1 12Hz-24Hz,coh F7-FC1 25Hz-29Hz,coh F7-FC2 1Hz-3Hz,coh F7-FC2 4Hz-7Hz,coh F7-FC2 8Hz-11Hz,coh F7-FC2 12Hz-24Hz,coh F7-FC2 25Hz-29Hz,coh F7-FC6 1Hz-3Hz,coh F7-FC6 4Hz-7Hz,coh F7-FC6 8Hz-11Hz,coh F7-FC6 12Hz-24Hz,coh F7-FC6 25Hz-29Hz,coh F7-T7 1Hz-3Hz,coh F7-T7 4Hz-7Hz,coh F7-T7 8Hz-11Hz,coh F7-T7 12Hz-24Hz,coh F7-T7 25Hz-29Hz,coh F7-C3 1Hz-3Hz,coh F7-C3 4Hz-7Hz,coh F7-C3 8Hz-11Hz,coh F7-C3 12Hz-24Hz,coh F7-C3 25Hz-29Hz,coh F7-Cz 1Hz-3Hz,coh F7-Cz 4Hz-7Hz,coh F7-Cz 8Hz-11Hz,coh F7-Cz 12Hz-24Hz,coh F7-Cz 25Hz-29Hz,coh F7-C4 1Hz-3Hz,coh F7-C4 4Hz-7Hz,coh F7-C4 8Hz-11Hz,coh F7-C4 12Hz-24Hz,coh F7-C4 25Hz-29Hz,coh F7-T8 1Hz-3Hz,coh F7-T8 4Hz-7Hz,coh F7-T8 8Hz-11Hz,coh F7-T8 12Hz-24Hz,coh F7-T8 25Hz-29Hz,coh F7-CP5 1Hz-3Hz,coh F7-CP5 4Hz-7Hz,coh F7-CP5 8Hz-11Hz,coh F7-CP5 12Hz-24Hz,coh F7-CP5 25Hz-29Hz,coh F7-CP1 1Hz-3Hz,coh F7-CP1 4Hz-7Hz,coh F7-CP1 8Hz-11Hz,coh F7-CP1 12Hz-24Hz,coh F7-CP1 25Hz-29Hz,coh F7-CP2 1Hz-3Hz,coh F7-CP2 4Hz-7Hz,coh F7-CP2 8Hz-11Hz,coh F7-CP2 12Hz-24Hz,coh F7-CP2 25Hz-29Hz,coh F7-CP6 1Hz-3Hz,coh F7-CP6 4Hz-7Hz,coh F7-CP6 8Hz-11Hz,coh F7-CP6 12Hz-24Hz,coh F7-CP6 25Hz-29Hz,coh F7-P7 1Hz-3Hz,coh F7-P7 4Hz-7Hz,coh F7-P7 8Hz-11Hz,coh F7-P7 12Hz-24Hz,coh F7-P7 25Hz-29Hz,coh F7-P3 1Hz-3Hz,coh F7-P3 4Hz-7Hz,coh F7-P3 8Hz-11Hz,coh F7-P3 12Hz-24Hz,coh F7-P3 25Hz-29Hz,coh F7-Pz 1Hz-3Hz,coh F7-Pz 4Hz-7Hz,coh F7-Pz 8Hz-11Hz,coh F7-Pz 12Hz-24Hz,coh F7-Pz 25Hz-29Hz,coh F7-P4 1Hz-3Hz,coh F7-P4 4Hz-7Hz,coh F7-P4 8Hz-11Hz,coh F7-P4 12Hz-24Hz,coh F7-P4 25Hz-29Hz,coh F7-P8 1Hz-3Hz,coh F7-P8 4Hz-7Hz,coh F7-P8 8Hz-11Hz,coh F7-P8 12Hz-24Hz,coh F7-P8 25Hz-29Hz,coh F7-POz 1Hz-3Hz,coh F7-POz 4Hz-7Hz,coh F7-POz 8Hz-11Hz,coh F7-POz 12Hz-24Hz,coh F7-POz 25Hz-29Hz,coh F7-O1 1Hz-3Hz,coh F7-O1 4Hz-7Hz,coh F7-O1 8Hz-11Hz,coh F7-O1 12Hz-24Hz,coh F7-O1 25Hz-29Hz,coh F7-Oz 1Hz-3Hz,coh F7-Oz 4Hz-7Hz,coh F7-Oz 8Hz-11Hz,coh F7-Oz 12Hz-24Hz,coh F7-Oz 25Hz-29Hz,coh F7-O2 1Hz-3Hz,coh F7-O2 4Hz-7Hz,coh F7-O2 8Hz-11Hz,coh F7-O2 12Hz-24Hz,coh F7-O2 25Hz-29Hz,coh F3-Fz 1Hz-3Hz,coh F3-Fz 4Hz-7Hz,coh F3-Fz 8Hz-11Hz,coh F3-Fz 12Hz-24Hz,coh F3-Fz 25Hz-29Hz,coh F3-F4 1Hz-3Hz,coh F3-F4 4Hz-7Hz,coh F3-F4 8Hz-11Hz,coh F3-F4 12Hz-24Hz,coh F3-F4 25Hz-29Hz,coh F3-F8 1Hz-3Hz,coh F3-F8 4Hz-7Hz,coh F3-F8 8Hz-11Hz,coh F3-F8 12Hz-24Hz,coh F3-F8 25Hz-29Hz,coh F3-FC5 1Hz-3Hz,coh F3-FC5 4Hz-7Hz,coh F3-FC5 8Hz-11Hz,coh F3-FC5 12Hz-24Hz,coh F3-FC5 25Hz-29Hz,coh F3-FC1 1Hz-3Hz,coh F3-FC1 4Hz-7Hz,coh F3-FC1 8Hz-11Hz,coh F3-FC1 12Hz-24Hz,coh F3-FC1 25Hz-29Hz,coh F3-FC2 1Hz-3Hz,coh F3-FC2 4Hz-7Hz,coh F3-FC2 8Hz-11Hz,coh F3-FC2 12Hz-24Hz,coh F3-FC2 25Hz-29Hz,coh F3-FC6 1Hz-3Hz,coh F3-FC6 4Hz-7Hz,coh F3-FC6 8Hz-11Hz,coh F3-FC6 12Hz-24Hz,coh F3-FC6 25Hz-29Hz,coh F3-T7 1Hz-3Hz,coh F3-T7 4Hz-7Hz,coh F3-T7 8Hz-11Hz,coh F3-T7 12Hz-24Hz,coh F3-T7 25Hz-29Hz,coh F3-C3 1Hz-3Hz,coh F3-C3 4Hz-7Hz,coh F3-C3 8Hz-11Hz,coh F3-C3 12Hz-24Hz,coh F3-C3 25Hz-29Hz,coh F3-Cz 1Hz-3Hz,coh F3-Cz 4Hz-7Hz,coh F3-Cz 8Hz-11Hz,coh F3-Cz 12Hz-24Hz,coh F3-Cz 25Hz-29Hz,coh F3-C4 1Hz-3Hz,coh F3-C4 4Hz-7Hz,coh F3-C4 8Hz-11Hz,coh F3-C4 12Hz-24Hz,coh F3-C4 25Hz-29Hz,coh F3-T8 1Hz-3Hz,coh F3-T8 4Hz-7Hz,coh F3-T8 8Hz-11Hz,coh F3-T8 12Hz-24Hz,coh F3-T8 25Hz-29Hz,coh F3-CP5 1Hz-3Hz,coh F3-CP5 4Hz-7Hz,coh F3-CP5 8Hz-11Hz,coh F3-CP5 12Hz-24Hz,coh F3-CP5 25Hz-29Hz,coh F3-CP1 1Hz-3Hz,coh F3-CP1 4Hz-7Hz,coh F3-CP1 8Hz-11Hz,coh F3-CP1 12Hz-24Hz,coh F3-CP1 25Hz-29Hz,coh F3-CP2 1Hz-3Hz,coh F3-CP2 4Hz-7Hz,coh F3-CP2 8Hz-11Hz,coh F3-CP2 12Hz-24Hz,coh F3-CP2 25Hz-29Hz,coh F3-CP6 1Hz-3Hz,coh F3-CP6 4Hz-7Hz,coh F3-CP6 8Hz-11Hz,coh F3-CP6 12Hz-24Hz,coh F3-CP6 25Hz-29Hz,coh F3-P7 1Hz-3Hz,coh F3-P7 4Hz-7Hz,coh F3-P7 8Hz-11Hz,coh F3-P7 12Hz-24Hz,coh F3-P7 25Hz-29Hz,coh F3-P3 1Hz-3Hz,coh F3-P3 4Hz-7Hz,coh F3-P3 8Hz-11Hz,coh F3-P3 12Hz-24Hz,coh F3-P3 25Hz-29Hz,coh F3-Pz 1Hz-3Hz,coh F3-Pz 4Hz-7Hz,coh F3-Pz 8Hz-11Hz,coh F3-Pz 12Hz-24Hz,coh F3-Pz 25Hz-29Hz,coh F3-P4 1Hz-3Hz,coh F3-P4 4Hz-7Hz,coh F3-P4 8Hz-11Hz,coh F3-P4 12Hz-24Hz,coh F3-P4 25Hz-29Hz,coh F3-P8 1Hz-3Hz,coh F3-P8 4Hz-7Hz,coh F3-P8 8Hz-11Hz,coh F3-P8 12Hz-24Hz,coh F3-P8 25Hz-29Hz,coh F3-POz 1Hz-3Hz,coh F3-POz 4Hz-7Hz,coh F3-POz 8Hz-11Hz,coh F3-POz 12Hz-24Hz,coh F3-POz 25Hz-29Hz,coh F3-O1 1Hz-3Hz,coh F3-O1 4Hz-7Hz,coh F3-O1 8Hz-11Hz,coh F3-O1 12Hz-24Hz,coh F3-O1 25Hz-29Hz,coh F3-Oz 1Hz-3Hz,coh F3-Oz 4Hz-7Hz,coh F3-Oz 8Hz-11Hz,coh F3-Oz 12Hz-24Hz,coh F3-Oz 25Hz-29Hz,coh F3-O2 1Hz-3Hz,coh F3-O2 4Hz-7Hz,coh F3-O2 8Hz-11Hz,coh F3-O2 12Hz-24Hz,coh F3-O2 25Hz-29Hz,coh Fz-F4 1Hz-3Hz,coh Fz-F4 4Hz-7Hz,coh Fz-F4 8Hz-11Hz,coh Fz-F4 12Hz-24Hz,coh Fz-F4 25Hz-29Hz,coh Fz-F8 1Hz-3Hz,coh Fz-F8 4Hz-7Hz,coh Fz-F8 8Hz-11Hz,coh Fz-F8 12Hz-24Hz,coh Fz-F8 25Hz-29Hz,coh Fz-FC5 1Hz-3Hz,coh Fz-FC5 4Hz-7Hz,coh Fz-FC5 8Hz-11Hz,coh Fz-FC5 12Hz-24Hz,coh Fz-FC5 25Hz-29Hz,coh Fz-FC1 1Hz-3Hz,coh Fz-FC1 4Hz-7Hz,coh Fz-FC1 8Hz-11Hz,coh Fz-FC1 12Hz-24Hz,coh Fz-FC1 25Hz-29Hz,coh Fz-FC2 1Hz-3Hz,coh Fz-FC2 4Hz-7Hz,coh Fz-FC2 8Hz-11Hz,coh Fz-FC2 12Hz-24Hz,coh Fz-FC2 25Hz-29Hz,coh Fz-FC6 1Hz-3Hz,coh Fz-FC6 4Hz-7Hz,coh Fz-FC6 8Hz-11Hz,coh Fz-FC6 12Hz-24Hz,coh Fz-FC6 25Hz-29Hz,coh Fz-T7 1Hz-3Hz,coh Fz-T7 4Hz-7Hz,coh Fz-T7 8Hz-11Hz,coh Fz-T7 12Hz-24Hz,coh Fz-T7 25Hz-29Hz,coh Fz-C3 1Hz-3Hz,coh Fz-C3 4Hz-7Hz,coh Fz-C3 8Hz-11Hz,coh Fz-C3 12Hz-24Hz,coh Fz-C3 25Hz-29Hz,coh Fz-Cz 1Hz-3Hz,coh Fz-Cz 4Hz-7Hz,coh Fz-Cz 8Hz-11Hz,coh Fz-Cz 12Hz-24Hz,coh Fz-Cz 25Hz-29Hz,coh Fz-C4 1Hz-3Hz,coh Fz-C4 4Hz-7Hz,coh Fz-C4 8Hz-11Hz,coh Fz-C4 12Hz-24Hz,coh Fz-C4 25Hz-29Hz,coh Fz-T8 1Hz-3Hz,coh Fz-T8 4Hz-7Hz,coh Fz-T8 8Hz-11Hz,coh Fz-T8 12Hz-24Hz,coh Fz-T8 25Hz-29Hz,coh Fz-CP5 1Hz-3Hz,coh Fz-CP5 4Hz-7Hz,coh Fz-CP5 8Hz-11Hz,coh Fz-CP5 12Hz-24Hz,coh Fz-CP5 25Hz-29Hz,coh Fz-CP1 1Hz-3Hz,coh Fz-CP1 4Hz-7Hz,coh Fz-CP1 8Hz-11Hz,coh Fz-CP1 12Hz-24Hz,coh Fz-CP1 25Hz-29Hz,coh Fz-CP2 1Hz-3Hz,coh Fz-CP2 4Hz-7Hz,coh Fz-CP2 8Hz-11Hz,coh Fz-CP2 12Hz-24Hz,coh Fz-CP2 25Hz-29Hz,coh Fz-CP6 1Hz-3Hz,coh Fz-CP6 4Hz-7Hz,coh Fz-CP6 8Hz-11Hz,coh Fz-CP6 12Hz-24Hz,coh Fz-CP6 25Hz-29Hz,coh Fz-P7 1Hz-3Hz,coh Fz-P7 4Hz-7Hz,coh Fz-P7 8Hz-11Hz,coh Fz-P7 12Hz-24Hz,coh Fz-P7 25Hz-29Hz,coh Fz-P3 1Hz-3Hz,coh Fz-P3 4Hz-7Hz,coh Fz-P3 8Hz-11Hz,coh Fz-P3 12Hz-24Hz,coh Fz-P3 25Hz-29Hz,coh Fz-Pz 1Hz-3Hz,coh Fz-Pz 4Hz-7Hz,coh Fz-Pz 8Hz-11Hz,coh Fz-Pz 12Hz-24Hz,coh Fz-Pz 25Hz-29Hz,coh Fz-P4 1Hz-3Hz,coh Fz-P4 4Hz-7Hz,coh Fz-P4 8Hz-11Hz,coh Fz-P4 12Hz-24Hz,coh Fz-P4 25Hz-29Hz,coh Fz-P8 1Hz-3Hz,coh Fz-P8 4Hz-7Hz,coh Fz-P8 8Hz-11Hz,coh Fz-P8 12Hz-24Hz,coh Fz-P8 25Hz-29Hz,coh Fz-POz 1Hz-3Hz,coh Fz-POz 4Hz-7Hz,coh Fz-POz 8Hz-11Hz,coh Fz-POz 12Hz-24Hz,coh Fz-POz 25Hz-29Hz,coh Fz-O1 1Hz-3Hz,coh Fz-O1 4Hz-7Hz,coh Fz-O1 8Hz-11Hz,coh Fz-O1 12Hz-24Hz,coh Fz-O1 25Hz-29Hz,coh Fz-Oz 1Hz-3Hz,coh Fz-Oz 4Hz-7Hz,coh Fz-Oz 8Hz-11Hz,coh Fz-Oz 12Hz-24Hz,coh Fz-Oz 25Hz-29Hz,coh Fz-O2 1Hz-3Hz,coh Fz-O2 4Hz-7Hz,coh Fz-O2 8Hz-11Hz,coh Fz-O2 12Hz-24Hz,coh Fz-O2 25Hz-29Hz,coh F4-F8 1Hz-3Hz,coh F4-F8 4Hz-7Hz,coh F4-F8 8Hz-11Hz,coh F4-F8 12Hz-24Hz,coh F4-F8 25Hz-29Hz,coh F4-FC5 1Hz-3Hz,coh F4-FC5 4Hz-7Hz,coh F4-FC5 8Hz-11Hz,coh F4-FC5 12Hz-24Hz,coh F4-FC5 25Hz-29Hz,coh F4-FC1 1Hz-3Hz,coh F4-FC1 4Hz-7Hz,coh F4-FC1 8Hz-11Hz,coh F4-FC1 12Hz-24Hz,coh F4-FC1 25Hz-29Hz,coh F4-FC2 1Hz-3Hz,coh F4-FC2 4Hz-7Hz,coh F4-FC2 8Hz-11Hz,coh F4-FC2 12Hz-24Hz,coh F4-FC2 25Hz-29Hz,coh F4-FC6 1Hz-3Hz,coh F4-FC6 4Hz-7Hz,coh F4-FC6 8Hz-11Hz,coh F4-FC6 12Hz-24Hz,coh F4-FC6 25Hz-29Hz,coh F4-T7 1Hz-3Hz,coh F4-T7 4Hz-7Hz,coh F4-T7 8Hz-11Hz,coh F4-T7 12Hz-24Hz,coh F4-T7 25Hz-29Hz,coh F4-C3 1Hz-3Hz,coh F4-C3 4Hz-7Hz,coh F4-C3 8Hz-11Hz,coh F4-C3 12Hz-24Hz,coh F4-C3 25Hz-29Hz,coh F4-Cz 1Hz-3Hz,coh F4-Cz 4Hz-7Hz,coh F4-Cz 8Hz-11Hz,coh F4-Cz 12Hz-24Hz,coh F4-Cz 25Hz-29Hz,coh F4-C4 1Hz-3Hz,coh F4-C4 4Hz-7Hz,coh F4-C4 8Hz-11Hz,coh F4-C4 12Hz-24Hz,coh F4-C4 25Hz-29Hz,coh F4-T8 1Hz-3Hz,coh F4-T8 4Hz-7Hz,coh F4-T8 8Hz-11Hz,coh F4-T8 12Hz-24Hz,coh F4-T8 25Hz-29Hz,coh F4-CP5 1Hz-3Hz,coh F4-CP5 4Hz-7Hz,coh F4-CP5 8Hz-11Hz,coh F4-CP5 12Hz-24Hz,coh F4-CP5 25Hz-29Hz,coh F4-CP1 1Hz-3Hz,coh F4-CP1 4Hz-7Hz,coh F4-CP1 8Hz-11Hz,coh F4-CP1 12Hz-24Hz,coh F4-CP1 25Hz-29Hz,coh F4-CP2 1Hz-3Hz,coh F4-CP2 4Hz-7Hz,coh F4-CP2 8Hz-11Hz,coh F4-CP2 12Hz-24Hz,coh F4-CP2 25Hz-29Hz,coh F4-CP6 1Hz-3Hz,coh F4-CP6 4Hz-7Hz,coh F4-CP6 8Hz-11Hz,coh F4-CP6 12Hz-24Hz,coh F4-CP6 25Hz-29Hz,coh F4-P7 1Hz-3Hz,coh F4-P7 4Hz-7Hz,coh F4-P7 8Hz-11Hz,coh F4-P7 12Hz-24Hz,coh F4-P7 25Hz-29Hz,coh F4-P3 1Hz-3Hz,coh F4-P3 4Hz-7Hz,coh F4-P3 8Hz-11Hz,coh F4-P3 12Hz-24Hz,coh F4-P3 25Hz-29Hz,coh F4-Pz 1Hz-3Hz,coh F4-Pz 4Hz-7Hz,coh F4-Pz 8Hz-11Hz,coh F4-Pz 12Hz-24Hz,coh F4-Pz 25Hz-29Hz,coh F4-P4 1Hz-3Hz,coh F4-P4 4Hz-7Hz,coh F4-P4 8Hz-11Hz,coh F4-P4 12Hz-24Hz,coh F4-P4 25Hz-29Hz,coh F4-P8 1Hz-3Hz,coh F4-P8 4Hz-7Hz,coh F4-P8 8Hz-11Hz,coh F4-P8 12Hz-24Hz,coh F4-P8 25Hz-29Hz,coh F4-POz 1Hz-3Hz,coh F4-POz 4Hz-7Hz,coh F4-POz 8Hz-11Hz,coh F4-POz 12Hz-24Hz,coh F4-POz 25Hz-29Hz,coh F4-O1 1Hz-3Hz,coh F4-O1 4Hz-7Hz,coh F4-O1 8Hz-11Hz,coh F4-O1 12Hz-24Hz,coh F4-O1 25Hz-29Hz,coh F4-Oz 1Hz-3Hz,coh F4-Oz 4Hz-7Hz,coh F4-Oz 8Hz-11Hz,coh F4-Oz 12Hz-24Hz,coh F4-Oz 25Hz-29Hz,coh F4-O2 1Hz-3Hz,coh F4-O2 4Hz-7Hz,coh F4-O2 8Hz-11Hz,coh F4-O2 12Hz-24Hz,coh F4-O2 25Hz-29Hz,coh F8-FC5 1Hz-3Hz,coh F8-FC5 4Hz-7Hz,coh F8-FC5 8Hz-11Hz,coh F8-FC5 12Hz-24Hz,coh F8-FC5 25Hz-29Hz,coh F8-FC1 1Hz-3Hz,coh F8-FC1 4Hz-7Hz,coh F8-FC1 8Hz-11Hz,coh F8-FC1 12Hz-24Hz,coh F8-FC1 25Hz-29Hz,coh F8-FC2 1Hz-3Hz,coh F8-FC2 4Hz-7Hz,coh F8-FC2 8Hz-11Hz,coh F8-FC2 12Hz-24Hz,coh F8-FC2 25Hz-29Hz,coh F8-FC6 1Hz-3Hz,coh F8-FC6 4Hz-7Hz,coh F8-FC6 8Hz-11Hz,coh F8-FC6 12Hz-24Hz,coh F8-FC6 25Hz-29Hz,coh F8-T7 1Hz-3Hz,coh F8-T7 4Hz-7Hz,coh F8-T7 8Hz-11Hz,coh F8-T7 12Hz-24Hz,coh F8-T7 25Hz-29Hz,coh F8-C3 1Hz-3Hz,coh F8-C3 4Hz-7Hz,coh F8-C3 8Hz-11Hz,coh F8-C3 12Hz-24Hz,coh F8-C3 25Hz-29Hz,coh F8-Cz 1Hz-3Hz,coh F8-Cz 4Hz-7Hz,coh F8-Cz 8Hz-11Hz,coh F8-Cz 12Hz-24Hz,coh F8-Cz 25Hz-29Hz,coh F8-C4 1Hz-3Hz,coh F8-C4 4Hz-7Hz,coh F8-C4 8Hz-11Hz,coh F8-C4 12Hz-24Hz,coh F8-C4 25Hz-29Hz,coh F8-T8 1Hz-3Hz,coh F8-T8 4Hz-7Hz,coh F8-T8 8Hz-11Hz,coh F8-T8 12Hz-24Hz,coh F8-T8 25Hz-29Hz,coh F8-CP5 1Hz-3Hz,coh F8-CP5 4Hz-7Hz,coh F8-CP5 8Hz-11Hz,coh F8-CP5 12Hz-24Hz,coh F8-CP5 25Hz-29Hz,coh F8-CP1 1Hz-3Hz,coh F8-CP1 4Hz-7Hz,coh F8-CP1 8Hz-11Hz,coh F8-CP1 12Hz-24Hz,coh F8-CP1 25Hz-29Hz,coh F8-CP2 1Hz-3Hz,coh F8-CP2 4Hz-7Hz,coh F8-CP2 8Hz-11Hz,coh F8-CP2 12Hz-24Hz,coh F8-CP2 25Hz-29Hz,coh F8-CP6 1Hz-3Hz,coh F8-CP6 4Hz-7Hz,coh F8-CP6 8Hz-11Hz,coh F8-CP6 12Hz-24Hz,coh F8-CP6 25Hz-29Hz,coh F8-P7 1Hz-3Hz,coh F8-P7 4Hz-7Hz,coh F8-P7 8Hz-11Hz,coh F8-P7 12Hz-24Hz,coh F8-P7 25Hz-29Hz,coh F8-P3 1Hz-3Hz,coh F8-P3 4Hz-7Hz,coh F8-P3 8Hz-11Hz,coh F8-P3 12Hz-24Hz,coh F8-P3 25Hz-29Hz,coh F8-Pz 1Hz-3Hz,coh F8-Pz 4Hz-7Hz,coh F8-Pz 8Hz-11Hz,coh F8-Pz 12Hz-24Hz,coh F8-Pz 25Hz-29Hz,coh F8-P4 1Hz-3Hz,coh F8-P4 4Hz-7Hz,coh F8-P4 8Hz-11Hz,coh F8-P4 12Hz-24Hz,coh F8-P4 25Hz-29Hz,coh F8-P8 1Hz-3Hz,coh F8-P8 4Hz-7Hz,coh F8-P8 8Hz-11Hz,coh F8-P8 12Hz-24Hz,coh F8-P8 25Hz-29Hz,coh F8-POz 1Hz-3Hz,coh F8-POz 4Hz-7Hz,coh F8-POz 8Hz-11Hz,coh F8-POz 12Hz-24Hz,coh F8-POz 25Hz-29Hz,coh F8-O1 1Hz-3Hz,coh F8-O1 4Hz-7Hz,coh F8-O1 8Hz-11Hz,coh F8-O1 12Hz-24Hz,coh F8-O1 25Hz-29Hz,coh F8-Oz 1Hz-3Hz,coh F8-Oz 4Hz-7Hz,coh F8-Oz 8Hz-11Hz,coh F8-Oz 12Hz-24Hz,coh F8-Oz 25Hz-29Hz,coh F8-O2 1Hz-3Hz,coh F8-O2 4Hz-7Hz,coh F8-O2 8Hz-11Hz,coh F8-O2 12Hz-24Hz,coh F8-O2 25Hz-29Hz,coh FC5-FC1 1Hz-3Hz,coh FC5-FC1 4Hz-7Hz,coh FC5-FC1 8Hz-11Hz,coh FC5-FC1 12Hz-24Hz,coh FC5-FC1 25Hz-29Hz,coh FC5-FC2 1Hz-3Hz,coh FC5-FC2 4Hz-7Hz,coh FC5-FC2 8Hz-11Hz,coh FC5-FC2 12Hz-24Hz,coh FC5-FC2 25Hz-29Hz,coh FC5-FC6 1Hz-3Hz,coh FC5-FC6 4Hz-7Hz,coh FC5-FC6 8Hz-11Hz,coh FC5-FC6 12Hz-24Hz,coh FC5-FC6 25Hz-29Hz,coh FC5-T7 1Hz-3Hz,coh FC5-T7 4Hz-7Hz,coh FC5-T7 8Hz-11Hz,coh FC5-T7 12Hz-24Hz,coh FC5-T7 25Hz-29Hz,coh FC5-C3 1Hz-3Hz,coh FC5-C3 4Hz-7Hz,coh FC5-C3 8Hz-11Hz,coh FC5-C3 12Hz-24Hz,coh FC5-C3 25Hz-29Hz,coh FC5-Cz 1Hz-3Hz,coh FC5-Cz 4Hz-7Hz,coh FC5-Cz 8Hz-11Hz,coh FC5-Cz 12Hz-24Hz,coh FC5-Cz 25Hz-29Hz,coh FC5-C4 1Hz-3Hz,coh FC5-C4 4Hz-7Hz,coh FC5-C4 8Hz-11Hz,coh FC5-C4 12Hz-24Hz,coh FC5-C4 25Hz-29Hz,coh FC5-T8 1Hz-3Hz,coh FC5-T8 4Hz-7Hz,coh FC5-T8 8Hz-11Hz,coh FC5-T8 12Hz-24Hz,coh FC5-T8 25Hz-29Hz,coh FC5-CP5 1Hz-3Hz,coh FC5-CP5 4Hz-7Hz,coh FC5-CP5 8Hz-11Hz,coh FC5-CP5 12Hz-24Hz,coh FC5-CP5 25Hz-29Hz,coh FC5-CP1 1Hz-3Hz,coh FC5-CP1 4Hz-7Hz,coh FC5-CP1 8Hz-11Hz,coh FC5-CP1 12Hz-24Hz,coh FC5-CP1 25Hz-29Hz,coh FC5-CP2 1Hz-3Hz,coh FC5-CP2 4Hz-7Hz,coh FC5-CP2 8Hz-11Hz,coh FC5-CP2 12Hz-24Hz,coh FC5-CP2 25Hz-29Hz,coh FC5-CP6 1Hz-3Hz,coh FC5-CP6 4Hz-7Hz,coh FC5-CP6 8Hz-11Hz,coh FC5-CP6 12Hz-24Hz,coh FC5-CP6 25Hz-29Hz,coh FC5-P7 1Hz-3Hz,coh FC5-P7 4Hz-7Hz,coh FC5-P7 8Hz-11Hz,coh FC5-P7 12Hz-24Hz,coh FC5-P7 25Hz-29Hz,coh FC5-P3 1Hz-3Hz,coh FC5-P3 4Hz-7Hz,coh FC5-P3 8Hz-11Hz,coh FC5-P3 12Hz-24Hz,coh FC5-P3 25Hz-29Hz,coh FC5-Pz 1Hz-3Hz,coh FC5-Pz 4Hz-7Hz,coh FC5-Pz 8Hz-11Hz,coh FC5-Pz 12Hz-24Hz,coh FC5-Pz 25Hz-29Hz,coh FC5-P4 1Hz-3Hz,coh FC5-P4 4Hz-7Hz,coh FC5-P4 8Hz-11Hz,coh FC5-P4 12Hz-24Hz,coh FC5-P4 25Hz-29Hz,coh FC5-P8 1Hz-3Hz,coh FC5-P8 4Hz-7Hz,coh FC5-P8 8Hz-11Hz,coh FC5-P8 12Hz-24Hz,coh FC5-P8 25Hz-29Hz,coh FC5-POz 1Hz-3Hz,coh FC5-POz 4Hz-7Hz,coh FC5-POz 8Hz-11Hz,coh FC5-POz 12Hz-24Hz,coh FC5-POz 25Hz-29Hz,coh FC5-O1 1Hz-3Hz,coh FC5-O1 4Hz-7Hz,coh FC5-O1 8Hz-11Hz,coh FC5-O1 12Hz-24Hz,coh FC5-O1 25Hz-29Hz,coh FC5-Oz 1Hz-3Hz,coh FC5-Oz 4Hz-7Hz,coh FC5-Oz 8Hz-11Hz,coh FC5-Oz 12Hz-24Hz,coh FC5-Oz 25Hz-29Hz,coh FC5-O2 1Hz-3Hz,coh FC5-O2 4Hz-7Hz,coh FC5-O2 8Hz-11Hz,coh FC5-O2 12Hz-24Hz,coh FC5-O2 25Hz-29Hz,coh FC1-FC2 1Hz-3Hz,coh FC1-FC2 4Hz-7Hz,coh FC1-FC2 8Hz-11Hz,coh FC1-FC2 12Hz-24Hz,coh FC1-FC2 25Hz-29Hz,coh FC1-FC6 1Hz-3Hz,coh FC1-FC6 4Hz-7Hz,coh FC1-FC6 8Hz-11Hz,coh FC1-FC6 12Hz-24Hz,coh FC1-FC6 25Hz-29Hz,coh FC1-T7 1Hz-3Hz,coh FC1-T7 4Hz-7Hz,coh FC1-T7 8Hz-11Hz,coh FC1-T7 12Hz-24Hz,coh FC1-T7 25Hz-29Hz,coh FC1-C3 1Hz-3Hz,coh FC1-C3 4Hz-7Hz,coh FC1-C3 8Hz-11Hz,coh FC1-C3 12Hz-24Hz,coh FC1-C3 25Hz-29Hz,coh FC1-Cz 1Hz-3Hz,coh FC1-Cz 4Hz-7Hz,coh FC1-Cz 8Hz-11Hz,coh FC1-Cz 12Hz-24Hz,coh FC1-Cz 25Hz-29Hz,coh FC1-C4 1Hz-3Hz,coh FC1-C4 4Hz-7Hz,coh FC1-C4 8Hz-11Hz,coh FC1-C4 12Hz-24Hz,coh FC1-C4 25Hz-29Hz,coh FC1-T8 1Hz-3Hz,coh FC1-T8 4Hz-7Hz,coh FC1-T8 8Hz-11Hz,coh FC1-T8 12Hz-24Hz,coh FC1-T8 25Hz-29Hz,coh FC1-CP5 1Hz-3Hz,coh FC1-CP5 4Hz-7Hz,coh FC1-CP5 8Hz-11Hz,coh FC1-CP5 12Hz-24Hz,coh FC1-CP5 25Hz-29Hz,coh FC1-CP1 1Hz-3Hz,coh FC1-CP1 4Hz-7Hz,coh FC1-CP1 8Hz-11Hz,coh FC1-CP1 12Hz-24Hz,coh FC1-CP1 25Hz-29Hz,coh FC1-CP2 1Hz-3Hz,coh FC1-CP2 4Hz-7Hz,coh FC1-CP2 8Hz-11Hz,coh FC1-CP2 12Hz-24Hz,coh FC1-CP2 25Hz-29Hz,coh FC1-CP6 1Hz-3Hz,coh FC1-CP6 4Hz-7Hz,coh FC1-CP6 8Hz-11Hz,coh FC1-CP6 12Hz-24Hz,coh FC1-CP6 25Hz-29Hz,coh FC1-P7 1Hz-3Hz,coh FC1-P7 4Hz-7Hz,coh FC1-P7 8Hz-11Hz,coh FC1-P7 12Hz-24Hz,coh FC1-P7 25Hz-29Hz,coh FC1-P3 1Hz-3Hz,coh FC1-P3 4Hz-7Hz,coh FC1-P3 8Hz-11Hz,coh FC1-P3 12Hz-24Hz,coh FC1-P3 25Hz-29Hz,coh FC1-Pz 1Hz-3Hz,coh FC1-Pz 4Hz-7Hz,coh FC1-Pz 8Hz-11Hz,coh FC1-Pz 12Hz-24Hz,coh FC1-Pz 25Hz-29Hz,coh FC1-P4 1Hz-3Hz,coh FC1-P4 4Hz-7Hz,coh FC1-P4 8Hz-11Hz,coh FC1-P4 12Hz-24Hz,coh FC1-P4 25Hz-29Hz,coh FC1-P8 1Hz-3Hz,coh FC1-P8 4Hz-7Hz,coh FC1-P8 8Hz-11Hz,coh FC1-P8 12Hz-24Hz,coh FC1-P8 25Hz-29Hz,coh FC1-POz 1Hz-3Hz,coh FC1-POz 4Hz-7Hz,coh FC1-POz 8Hz-11Hz,coh FC1-POz 12Hz-24Hz,coh FC1-POz 25Hz-29Hz,coh FC1-O1 1Hz-3Hz,coh FC1-O1 4Hz-7Hz,coh FC1-O1 8Hz-11Hz,coh FC1-O1 12Hz-24Hz,coh FC1-O1 25Hz-29Hz,coh FC1-Oz 1Hz-3Hz,coh FC1-Oz 4Hz-7Hz,coh FC1-Oz 8Hz-11Hz,coh FC1-Oz 12Hz-24Hz,coh FC1-Oz 25Hz-29Hz,coh FC1-O2 1Hz-3Hz,coh FC1-O2 4Hz-7Hz,coh FC1-O2 8Hz-11Hz,coh FC1-O2 12Hz-24Hz,coh FC1-O2 25Hz-29Hz,coh FC2-FC6 1Hz-3Hz,coh FC2-FC6 4Hz-7Hz,coh FC2-FC6 8Hz-11Hz,coh FC2-FC6 12Hz-24Hz,coh FC2-FC6 25Hz-29Hz,coh FC2-T7 1Hz-3Hz,coh FC2-T7 4Hz-7Hz,coh FC2-T7 8Hz-11Hz,coh FC2-T7 12Hz-24Hz,coh FC2-T7 25Hz-29Hz,coh FC2-C3 1Hz-3Hz,coh FC2-C3 4Hz-7Hz,coh FC2-C3 8Hz-11Hz,coh FC2-C3 12Hz-24Hz,coh FC2-C3 25Hz-29Hz,coh FC2-Cz 1Hz-3Hz,coh FC2-Cz 4Hz-7Hz,coh FC2-Cz 8Hz-11Hz,coh FC2-Cz 12Hz-24Hz,coh FC2-Cz 25Hz-29Hz,coh FC2-C4 1Hz-3Hz,coh FC2-C4 4Hz-7Hz,coh FC2-C4 8Hz-11Hz,coh FC2-C4 12Hz-24Hz,coh FC2-C4 25Hz-29Hz,coh FC2-T8 1Hz-3Hz,coh FC2-T8 4Hz-7Hz,coh FC2-T8 8Hz-11Hz,coh FC2-T8 12Hz-24Hz,coh FC2-T8 25Hz-29Hz,coh FC2-CP5 1Hz-3Hz,coh FC2-CP5 4Hz-7Hz,coh FC2-CP5 8Hz-11Hz,coh FC2-CP5 12Hz-24Hz,coh FC2-CP5 25Hz-29Hz,coh FC2-CP1 1Hz-3Hz,coh FC2-CP1 4Hz-7Hz,coh FC2-CP1 8Hz-11Hz,coh FC2-CP1 12Hz-24Hz,coh FC2-CP1 25Hz-29Hz,coh FC2-CP2 1Hz-3Hz,coh FC2-CP2 4Hz-7Hz,coh FC2-CP2 8Hz-11Hz,coh FC2-CP2 12Hz-24Hz,coh FC2-CP2 25Hz-29Hz,coh FC2-CP6 1Hz-3Hz,coh FC2-CP6 4Hz-7Hz,coh FC2-CP6 8Hz-11Hz,coh FC2-CP6 12Hz-24Hz,coh FC2-CP6 25Hz-29Hz,coh FC2-P7 1Hz-3Hz,coh FC2-P7 4Hz-7Hz,coh FC2-P7 8Hz-11Hz,coh FC2-P7 12Hz-24Hz,coh FC2-P7 25Hz-29Hz,coh FC2-P3 1Hz-3Hz,coh FC2-P3 4Hz-7Hz,coh FC2-P3 8Hz-11Hz,coh FC2-P3 12Hz-24Hz,coh FC2-P3 25Hz-29Hz,coh FC2-Pz 1Hz-3Hz,coh FC2-Pz 4Hz-7Hz,coh FC2-Pz 8Hz-11Hz,coh FC2-Pz 12Hz-24Hz,coh FC2-Pz 25Hz-29Hz,coh FC2-P4 1Hz-3Hz,coh FC2-P4 4Hz-7Hz,coh FC2-P4 8Hz-11Hz,coh FC2-P4 12Hz-24Hz,coh FC2-P4 25Hz-29Hz,coh FC2-P8 1Hz-3Hz,coh FC2-P8 4Hz-7Hz,coh FC2-P8 8Hz-11Hz,coh FC2-P8 12Hz-24Hz,coh FC2-P8 25Hz-29Hz,coh FC2-POz 1Hz-3Hz,coh FC2-POz 4Hz-7Hz,coh FC2-POz 8Hz-11Hz,coh FC2-POz 12Hz-24Hz,coh FC2-POz 25Hz-29Hz,coh FC2-O1 1Hz-3Hz,coh FC2-O1 4Hz-7Hz,coh FC2-O1 8Hz-11Hz,coh FC2-O1 12Hz-24Hz,coh FC2-O1 25Hz-29Hz,coh FC2-Oz 1Hz-3Hz,coh FC2-Oz 4Hz-7Hz,coh FC2-Oz 8Hz-11Hz,coh FC2-Oz 12Hz-24Hz,coh FC2-Oz 25Hz-29Hz,coh FC2-O2 1Hz-3Hz,coh FC2-O2 4Hz-7Hz,coh FC2-O2 8Hz-11Hz,coh FC2-O2 12Hz-24Hz,coh FC2-O2 25Hz-29Hz,coh FC6-T7 1Hz-3Hz,coh FC6-T7 4Hz-7Hz,coh FC6-T7 8Hz-11Hz,coh FC6-T7 12Hz-24Hz,coh FC6-T7 25Hz-29Hz,coh FC6-C3 1Hz-3Hz,coh FC6-C3 4Hz-7Hz,coh FC6-C3 8Hz-11Hz,coh FC6-C3 12Hz-24Hz,coh FC6-C3 25Hz-29Hz,coh FC6-Cz 1Hz-3Hz,coh FC6-Cz 4Hz-7Hz,coh FC6-Cz 8Hz-11Hz,coh FC6-Cz 12Hz-24Hz,coh FC6-Cz 25Hz-29Hz,coh FC6-C4 1Hz-3Hz,coh FC6-C4 4Hz-7Hz,coh FC6-C4 8Hz-11Hz,coh FC6-C4 12Hz-24Hz,coh FC6-C4 25Hz-29Hz,coh FC6-T8 1Hz-3Hz,coh FC6-T8 4Hz-7Hz,coh FC6-T8 8Hz-11Hz,coh FC6-T8 12Hz-24Hz,coh FC6-T8 25Hz-29Hz,coh FC6-CP5 1Hz-3Hz,coh FC6-CP5 4Hz-7Hz,coh FC6-CP5 8Hz-11Hz,coh FC6-CP5 12Hz-24Hz,coh FC6-CP5 25Hz-29Hz,coh FC6-CP1 1Hz-3Hz,coh FC6-CP1 4Hz-7Hz,coh FC6-CP1 8Hz-11Hz,coh FC6-CP1 12Hz-24Hz,coh FC6-CP1 25Hz-29Hz,coh FC6-CP2 1Hz-3Hz,coh FC6-CP2 4Hz-7Hz,coh FC6-CP2 8Hz-11Hz,coh FC6-CP2 12Hz-24Hz,coh FC6-CP2 25Hz-29Hz,coh FC6-CP6 1Hz-3Hz,coh FC6-CP6 4Hz-7Hz,coh FC6-CP6 8Hz-11Hz,coh FC6-CP6 12Hz-24Hz,coh FC6-CP6 25Hz-29Hz,coh FC6-P7 1Hz-3Hz,coh FC6-P7 4Hz-7Hz,coh FC6-P7 8Hz-11Hz,coh FC6-P7 12Hz-24Hz,coh FC6-P7 25Hz-29Hz,coh FC6-P3 1Hz-3Hz,coh FC6-P3 4Hz-7Hz,coh FC6-P3 8Hz-11Hz,coh FC6-P3 12Hz-24Hz,coh FC6-P3 25Hz-29Hz,coh FC6-Pz 1Hz-3Hz,coh FC6-Pz 4Hz-7Hz,coh FC6-Pz 8Hz-11Hz,coh FC6-Pz 12Hz-24Hz,coh FC6-Pz 25Hz-29Hz,coh FC6-P4 1Hz-3Hz,coh FC6-P4 4Hz-7Hz,coh FC6-P4 8Hz-11Hz,coh FC6-P4 12Hz-24Hz,coh FC6-P4 25Hz-29Hz,coh FC6-P8 1Hz-3Hz,coh FC6-P8 4Hz-7Hz,coh FC6-P8 8Hz-11Hz,coh FC6-P8 12Hz-24Hz,coh FC6-P8 25Hz-29Hz,coh FC6-POz 1Hz-3Hz,coh FC6-POz 4Hz-7Hz,coh FC6-POz 8Hz-11Hz,coh FC6-POz 12Hz-24Hz,coh FC6-POz 25Hz-29Hz,coh FC6-O1 1Hz-3Hz,coh FC6-O1 4Hz-7Hz,coh FC6-O1 8Hz-11Hz,coh FC6-O1 12Hz-24Hz,coh FC6-O1 25Hz-29Hz,coh FC6-Oz 1Hz-3Hz,coh FC6-Oz 4Hz-7Hz,coh FC6-Oz 8Hz-11Hz,coh FC6-Oz 12Hz-24Hz,coh FC6-Oz 25Hz-29Hz,coh FC6-O2 1Hz-3Hz,coh FC6-O2 4Hz-7Hz,coh FC6-O2 8Hz-11Hz,coh FC6-O2 12Hz-24Hz,coh FC6-O2 25Hz-29Hz,coh T7-C3 1Hz-3Hz,coh T7-C3 4Hz-7Hz,coh T7-C3 8Hz-11Hz,coh T7-C3 12Hz-24Hz,coh T7-C3 25Hz-29Hz,coh T7-Cz 1Hz-3Hz,coh T7-Cz 4Hz-7Hz,coh T7-Cz 8Hz-11Hz,coh T7-Cz 12Hz-24Hz,coh T7-Cz 25Hz-29Hz,coh T7-C4 1Hz-3Hz,coh T7-C4 4Hz-7Hz,coh T7-C4 8Hz-11Hz,coh T7-C4 12Hz-24Hz,coh T7-C4 25Hz-29Hz,coh T7-T8 1Hz-3Hz,coh T7-T8 4Hz-7Hz,coh T7-T8 8Hz-11Hz,coh T7-T8 12Hz-24Hz,coh T7-T8 25Hz-29Hz,coh T7-CP5 1Hz-3Hz,coh T7-CP5 4Hz-7Hz,coh T7-CP5 8Hz-11Hz,coh T7-CP5 12Hz-24Hz,coh T7-CP5 25Hz-29Hz,coh T7-CP1 1Hz-3Hz,coh T7-CP1 4Hz-7Hz,coh T7-CP1 8Hz-11Hz,coh T7-CP1 12Hz-24Hz,coh T7-CP1 25Hz-29Hz,coh T7-CP2 1Hz-3Hz,coh T7-CP2 4Hz-7Hz,coh T7-CP2 8Hz-11Hz,coh T7-CP2 12Hz-24Hz,coh T7-CP2 25Hz-29Hz,coh T7-CP6 1Hz-3Hz,coh T7-CP6 4Hz-7Hz,coh T7-CP6 8Hz-11Hz,coh T7-CP6 12Hz-24Hz,coh T7-CP6 25Hz-29Hz,coh T7-P7 1Hz-3Hz,coh T7-P7 4Hz-7Hz,coh T7-P7 8Hz-11Hz,coh T7-P7 12Hz-24Hz,coh T7-P7 25Hz-29Hz,coh T7-P3 1Hz-3Hz,coh T7-P3 4Hz-7Hz,coh T7-P3 8Hz-11Hz,coh T7-P3 12Hz-24Hz,coh T7-P3 25Hz-29Hz,coh T7-Pz 1Hz-3Hz,coh T7-Pz 4Hz-7Hz,coh T7-Pz 8Hz-11Hz,coh T7-Pz 12Hz-24Hz,coh T7-Pz 25Hz-29Hz,coh T7-P4 1Hz-3Hz,coh T7-P4 4Hz-7Hz,coh T7-P4 8Hz-11Hz,coh T7-P4 12Hz-24Hz,coh T7-P4 25Hz-29Hz,coh T7-P8 1Hz-3Hz,coh T7-P8 4Hz-7Hz,coh T7-P8 8Hz-11Hz,coh T7-P8 12Hz-24Hz,coh T7-P8 25Hz-29Hz,coh T7-POz 1Hz-3Hz,coh T7-POz 4Hz-7Hz,coh T7-POz 8Hz-11Hz,coh T7-POz 12Hz-24Hz,coh T7-POz 25Hz-29Hz,coh T7-O1 1Hz-3Hz,coh T7-O1 4Hz-7Hz,coh T7-O1 8Hz-11Hz,coh T7-O1 12Hz-24Hz,coh T7-O1 25Hz-29Hz,coh T7-Oz 1Hz-3Hz,coh T7-Oz 4Hz-7Hz,coh T7-Oz 8Hz-11Hz,coh T7-Oz 12Hz-24Hz,coh T7-Oz 25Hz-29Hz,coh T7-O2 1Hz-3Hz,coh T7-O2 4Hz-7Hz,coh T7-O2 8Hz-11Hz,coh T7-O2 12Hz-24Hz,coh T7-O2 25Hz-29Hz,coh C3-Cz 1Hz-3Hz,coh C3-Cz 4Hz-7Hz,coh C3-Cz 8Hz-11Hz,coh C3-Cz 12Hz-24Hz,coh C3-Cz 25Hz-29Hz,coh C3-C4 1Hz-3Hz,coh C3-C4 4Hz-7Hz,coh C3-C4 8Hz-11Hz,coh C3-C4 12Hz-24Hz,coh C3-C4 25Hz-29Hz,coh C3-T8 1Hz-3Hz,coh C3-T8 4Hz-7Hz,coh C3-T8 8Hz-11Hz,coh C3-T8 12Hz-24Hz,coh C3-T8 25Hz-29Hz,coh C3-CP5 1Hz-3Hz,coh C3-CP5 4Hz-7Hz,coh C3-CP5 8Hz-11Hz,coh C3-CP5 12Hz-24Hz,coh C3-CP5 25Hz-29Hz,coh C3-CP1 1Hz-3Hz,coh C3-CP1 4Hz-7Hz,coh C3-CP1 8Hz-11Hz... Output truncated.  Text exceeds maximum line length of 25,000 characters for Command Window display.

fclose(fileId);


% 
% 
% %initialize arrays
% queueLength = refreshRate * coherenceMemoryDuration;
% queueLengthRecip = 1/queueLength;
% powerArray = NaN(queueLength, freqCount, EEG.nbchan);
% powerSum = zeros(freqCount, EEG.nbchan);
% powerAvg = zeros(freqCount, EEG.nbchan);
% dotProdArray = NaN(queueLength, freqCount, channelPairCount);
% crossProdArray = NaN(queueLength, freqCount, channelPairCount);
% dotProdSum = zeros(freqCount, channelPairCount);
% crossProdSum = zeros(freqCount, channelPairCount);
% 
% coherenceIndex = 1;
% coherenceCount = 1;
% coherenceFull = false;
% fftStartIndex = 1;
% fftEndIndex = fftStartIndex + sampleRate*fftWindowDuration - 1;
% ffts = NaN(fftEndIndex, channelPairCount);
% meanFfts = NaN(length(lowFreq), EEG.nbchan);
% freqRange = 1:length(lowFreq);
% 
% fprintf('-');
% 
% 
% while(fftEndIndex < size(EEG.data,2))
%   fprintf('\b');
%   if(mod(coherenceCount, 4)==0)
%     fprintf('\\');
%   elseif(mod(coherenceCount, 4)==1)
%     fprintf('-');
%   elseif(mod(coherenceCount, 4)==2)
%     fprintf('/');
%   elseif(mod(coherenceCount, 4)==3)
%     fprintf('|');
%   end
%   if(mod(coherenceCount, 1000)==0)
%     percentDone =  100 * fftEndIndex / size(EEG.data,2);
%     fprintf('\n%s (%.3f%%)-', char(datetime), percentDone);
%   end
%   for chanCounter = 1:EEG.nbchan
% %     tic;
%     ffts(:, chanCounter) = fft(EEG.data(chanCounter, fftStartIndex:fftEndIndex));
% %     clock.fft = toc;
% %     tic;
%     for i = iRange
%       if(processFrequenciesInGroups)
%         myRange = (lowFreq(i)*fftWindowDuration+1):(highFreq(i)*fftWindowDuration+1);
%         meanFfts(i, chanCounter) = mean(ffts(myRange, chanCounter));
%       else
%         meanFfts(i, chanCounter) = ffts(i, chanCounter);
%         
%       end
%       if(coherenceFull)
%         powerSum(i, chanCounter) = powerSum(i, chanCounter) - powerArray(coherenceIndex, i, chanCounter);
%       end
%       x = real(meanFfts(i,chanCounter));
%       y = imag(meanFfts(i,chanCounter));
%       pow = 2*(x * x + y * y);
%       powerArray(coherenceIndex, i, chanCounter) = pow;
%       powerSum(i, chanCounter) = powerSum(i, chanCounter) + powerArray(coherenceIndex, i, chanCounter);
%       avgPow(i,chanCounter) = powerSum(i, chanCounter) * queueLengthRecip;
%     end
% %     clock.average = toc;
% 
%   end
%   
%   channelPairIndex = 1;
%   for chan1 = 1:EEG.nbchan
%     for chan2 = chan1+1:EEG.nbchan
%       if(length(chan2) > 0)
%         for i = iRange
% %           tic;
%           if(coherenceFull)
%             dotProdSum(i, channelPairIndex) = dotProdSum(i, channelPairIndex) ...
%               - dotProdArray(coherenceIndex, i, channelPairIndex);
%             crossProdSum(i, channelPairIndex) = crossProdSum(i, channelPairIndex) ...
%               - crossProdArray(coherenceIndex, i, channelPairIndex);
%           end
%           x1 = real(meanFfts(i,chan1));
%           y1 = imag(meanFfts(i,chan1));
%           x2 = real(meanFfts(i,chan2));
%           y2 = imag(meanFfts(i,chan2));
%           dot = 2*(x1 * x2 + y1 * y2);
%           cross = 2*(x1 * y2 - x2 * y1);
%           dotProdArray(coherenceIndex, i, channelPairIndex) = dot;
%           crossProdArray(coherenceIndex, i, channelPairIndex) = cross;
%           
%           dotProdSum(i, channelPairIndex) = dotProdSum(i, channelPairIndex)...
%             + dotProdArray(coherenceIndex, i, channelPairIndex);
%           crossProdSum(i, channelPairIndex) = crossProdSum(i, channelPairIndex)...
%             + crossProdArray(coherenceIndex, i, channelPairIndex);
%           if(coherenceFull)
%             avgPow1 = avgPow(i, chan1);
%             avgPow2 = avgPow(i, chan2);
%             avgDot = dotProdSum(i, channelPairIndex) * queueLengthRecip;
%             avgCross = crossProdSum(i, channelPairIndex) * queueLengthRecip;
%             channelPairs(channelPairIndex).coherence(coherenceCount, i) = ...
%               (avgDot*avgDot + avgCross*avgCross)...
%               / (avgPow1*avgPow2);
%             channels(chan1).absolutePower(coherenceCount, i) = log(avgPow1);
%             %need to fill in the last channel...
%             if(chan1 == 1)
%               channels(chan2).absolutePower(coherenceCount, i) = log(avgPow2);
%             end
% %           clock.coherence = toc;
% 
%           end
%         end
%         channelPairIndex = channelPairIndex + 1;
%       end
%     end
%     if(coherenceFull)
%       thisPowerSumRecip = 1 / sum(channels(chan1).absolutePower(coherenceCount, :));
%       channels(chan1).relativePower(coherenceCount, :) = channels(chan1).absolutePower(coherenceCount, :) .* thisPowerSumRecip;
%     end
%   end
%   %     %debug
%   %     if(coherenceFull)
%   %         avgPow1 = avgPow(i, chan1);
%   %     end
%   %     %end debug
%   coherenceCount = coherenceCount + 1;
%   coherenceIndex = coherenceIndex + 1;
%   if(coherenceIndex > queueLength)
%     coherenceIndex = 1;
%     coherenceFull = true;
%   end
%   fftStartIndex = fftStartIndex + (sampleRate / refreshRate);
%   fftEndIndex = fftStartIndex + sampleRate*fftWindowDuration - 1;
% end
% 
% x = ((1:size(channelPairs(1).coherence,1))-1) ./ refreshRate;
% remove = isnan(channelPairs(1).coherence(:,1));
% for i = 1:length(channelPairs)
%   channelPairs(i).coherence(remove,:) = [];
% end
% for i = 1:length(channels)
%   channels(i).absolutePower(remove,:)=[];
% end
% x(remove) = [];



toc;

