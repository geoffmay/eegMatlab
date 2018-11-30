function convolvedFourier = s5convolveFourierAtMarkers(fourierFolder, markerTimes)

files = dir(fourierFolder);
files = files(cellfun(@length, strfind({files.name}, '.mat')) > 0);
fileNumber = 1;

filePath = fullfile(fourierFolder, files(fileNumber).name);
data = load(filePath);

totalSamples = length(data.times) * length(files);


%sampleRate = 1 / (fourierEeg.x(2) - fourierEeg.x(1));
sampleRate = 1 / (data.times(2) - data.times(1));
%startTime = fourierEeg.x(1);
startTime = data.times(1);

% sampleCount = length(fourierEeg.x);
freqCount = length(data.freqInfo.lowFrequencies);
%signalCount = (length(fourierEeg.coh) + length(fourierEeg.channels)) * freqCount;
signalCount = (length(data.channelPairs) + length(data.channels)) * freqCount;

downsampleRate = 5;

convolvedFourier.signal = NaN(signalCount, length(markerTimes));
convolvedFourier.times = markerTimes;

signalNumber = 1;


for channelNumber = 1:length(data.channels)
  measureLabel = 'absolutePower';
  chanLabel = data.channels(channelNumber).label;
  for freqNumber = 1:freqCount
      %added
      fprintf('%s loading %d of %d\n', char(datetime), signalNumber, signalCount);
      fullSignal = NaN(totalSamples, 1);
      fullTimes = NaN(totalSamples, 1);
      sampleIndex = 1;
      for(fileNumber = 1:length(files ))
          data = load(fullfile(fourierFolder, files(fileNumber).name));
          piece = data.channels(channelNumber).absolutePower(:, freqNumber);
          endIndex = sampleIndex+length(piece) - 1;
          fullSignal(sampleIndex:endIndex) = piece;
          fullTimes(sampleIndex:endIndex) = data.times;
          sampleIndex = sampleIndex + length(piece);
      end
      remove = isnan(fullSignal);
      fullSignal(remove) = [];
      fullTimes(remove) = [];
      %end added
    freqLabel = sprintf('%d-%dHz', data.freqInfo.lowFrequencies(freqNumber), data.freqInfo.highFrequencies(freqNumber));
    label = sprintf('%s %s %s', measureLabel, chanLabel, freqLabel);
    fprintf('%s convolving %d of %d (%s)\n', char(datetime), signalNumber, signalCount, label);
%     sig = fourierEeg.channels(channelNumber).absolutePower(:, freqNumber);
%     cSig = s7convolveSignal(sig, sampleRate, downsampleRate, markerTimes, fourierEeg.x);
    cSig = s7convolveSignal(fullSignal, sampleRate, downsampleRate, markerTimes, fullTimes);
    fprintf('%s done\n', char(datetime));
    
    convolvedFourier.labels{signalNumber} = label;
    convolvedFourier.signals(signalNumber, :) = cSig;
    signalNumber = signalNumber + 1;
  end
end


for chanPairNumber = 1:length(data.channelPairs)
  measureLabel = 'coherence';
  chanLabel = data.channelPairs(chanPairNumber).label;
  for freqNumber = 1:freqCount
            %added
      fprintf('%s loading %d of %d\n', char(datetime), signalNumber, signalCount);
      fullSignal = NaN(totalSamples, 1);
      fullTimes = NaN(totalSamples, 1);
      sampleIndex = 1;
      for(fileNumber = 1:length(files ))
          data = load(fullfile(fourierFolder, files(fileNumber).name));
          piece = data.channelPairs(channelPairNumber).coherence(:, freqNumber);
          endIndex = sampleIndex+length(piece) - 1;
          fullSignal(sampleIndex:endIndex) = piece;
          fullTimes(sampleIndex:endIndex) = data.times;
          sampleIndex = sampleIndex + length(piece);
      end
      remove = isnan(fullSignal);
      fullSignal(remove) = [];
      fullTimes(remove) = [];
      %end added

    freqLabel = sprintf('%d-%dHz', data.freqInfo.lowFrequencies(freqNumber), data.freqInfo.highFrequencies(freqNumber));
    label = sprintf('%s %s %s', measureLabel, chanLabel, freqLabel);
    fprintf('%s convolving %d of %d (%s)\n', char(datetime), signalNumber, signalCount, label);
%     sig = fourierEeg.coh(chanPairNumber).coherence(:, freqNumber);
%     cSig = s7convolveSignal(sig, sampleRate, downsampleRate, markerTimes, fourierEeg.x);
    cSig = s7convolveSignal(fullSignal, sampleRate, downsampleRate, markerTimes, fullTimes);
    convolvedFourier.labels{signalNumber} = label;
    convolvedFourier.signals(signalNumber, :) = cSig;
    signalNumber = signalNumber + 1;
  end
end


