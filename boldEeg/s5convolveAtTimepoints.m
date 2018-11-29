function convolvedFourier = s5convolveFourierAtMarkers(fourierEeg, markerTimes)

%Signal: channel x timepoint matrix
%Markers: timepoints where we are interested


sampleRate = 1 / (fourierEeg.x(2) - fourierEeg.x(1));
startTime = fourierEeg.x(1);
sampleCount = length(fourierEeg.x);
freqCount = length(fourierEeg.freqInfo.lowFrequencies);
signalCount = (length(fourierEeg.coh) + length(fourierEeg.channels)) * freqCount;

downsampleRate = 5;

convolvedFourier.signal = NaN(signalCount, length(markerTimes));
convolvedFourier.times = markerTimes;

signalNumber = 1;


for channelNumber = 1:length(fourierEeg.channels)
  measureLabel = 'absolutePower';
  chanLabel = fourierEeg.channels(channelNumber).label;
  for freqNumber = 1:freqCount
    freqLabel = sprintf('%d-%dHz', fourierEeg.freqInfo.lowFrequencies(freqNumber), fourierEeg.freqInfo.highFrequencies(freqNumber));
    label = sprintf('%s %s %s', measureLabel, chanLabel, freqLabel);
    fprintf('%s convolving %d of %d (%s)\n', char(datetime), signalNumber, signalCount, label);
    sig = fourierEeg.channels(channelNumber).absolutePower(:, freqNumber);
    cSig = s7convolveSignal(sig, sampleRate, downsampleRate, markerTimes, fourierEeg.x);
    convolvedFourier.labels{signalNumber} = label;
    convolvedFourier.signals(signalNumber, :) = cSig;
    signalNumber = signalNumber + 1;
  end
end


for chanPairNumber = 1:length(fourierEeg.coh)
  measureLabel = 'coherence';
  chanLabel = fourierEeg.coh(chanPairNumber).label;
  for freqNumber = 1:freqCount
    freqLabel = sprintf('%d-%dHz', fourierEeg.freqInfo.lowFrequencies(freqNumber), fourierEeg.freqInfo.highFrequencies(freqNumber));
    label = sprintf('%s %s %s', measureLabel, chanLabel, freqLabel);
    fprintf('%s convolving %d of %d (%s)\n', char(datetime), signalNumber, signalCount, label);
    sig = fourierEeg.coh(chanPairNumber).coherence(:, freqNumber);
    cSig = s7convolveSignal(sig, sampleRate, downsampleRate, markerTimes, fourierEeg.x);
    convolvedFourier.labels{signalNumber} = label;
    convolvedFourier.signals(signalNumber, :) = cSig;
    signalNumber = signalNumber + 1;
  end
end


