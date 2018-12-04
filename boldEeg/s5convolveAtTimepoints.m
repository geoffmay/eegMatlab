function convolvedFourier = s5convolveFourierAtMarkers(fourierPath, markerTimes)

skipComputation = true;
tic;
chunkMaxRam = 1024 * 1024 * 1024 * 30;

files = dir(fourierPath);
files = files(cellfun(@length, strfind({files.name}, '.mat')) > 0);
fileNumber = 1;

filePath = fullfile(fourierPath, files(fileNumber).name);
testData = load(filePath);

totalSamples = length(testData.times) * length(files);

%sampleRate = 1 / (fourierEeg.x(2) - fourierEeg.x(1));
sampleRate = 1 / (testData.times(2) - testData.times(1));
%startTime = fourierEeg.x(1);
startTime = testData.times(1);

% sampleCount = length(fourierEeg.x);
freqCount = length(testData.freqInfo.lowFrequencies);
%signalCount = (length(fourierEeg.coh) + length(fourierEeg.channels)) * freqCount;
signalCount = (length(testData.channelPairs) + length(testData.channels)) * freqCount;


downsampleRate = 5;


chunkChannelCount = floor(chunkMaxRam / (totalSamples * 8));
% % don't downsample immediately, because of effects on chunk edges
%chunkChannelCount = floor(chunkMaxRam * downsampleRate / (totalSamples * 8));
%downsampledChunk = NaN(chunkChannelCount, ceil(totalSamples / downsampleRate));
if(~skipComputation)
chunk = NaN(chunkChannelCount, totalSamples);
end

%convolvedFourier.signal = NaN(signalCount, length(markerTimes));
convolvedFourier.times = markerTimes;

signalNumber = 1;


%for channelNumber = 1:length(data.channels)


%%power

channelNumber = 1;
chunkChannelIndex = 1;
oldChunkChannelIndex = 1;
while(channelNumber <= length(testData.channels))
  measureLabel = 'absolutePower';
  chanLabel = testData.channels(channelNumber).label;
  %for freqNumber = 1:freqCount
  freqNumber = 1;
  while(freqNumber <= freqCount && channelNumber <= length(testData.channels))
      %       fullSignal = NaN(totalSamples, 1);
      fullTimes = NaN(totalSamples, 1);
      %       sampleIndex = 1;
      oldChunkChannelIndex = chunkChannelIndex;
      oldFreqNumber = freqNumber;
      oldChanNumber = channelNumber;
      fprintf('%s convolving power, chunk %d-%d of %d\n', char(datetime), signalNumber, signalNumber + chunkChannelCount-1, signalCount);
      
      %to conserve ram, copy saved data to array in chunks
      for(fileNumber = 1:length(files ))

          %piece = data.channels(channelNumber).absolutePower(:, freqNumber);
          chunkChannelIndex = oldChunkChannelIndex;
          freqNumber = oldFreqNumber;
          channelNumber = oldChanNumber;
          fprintf('%s file %d of %d, chunk %d-%d\n', char(datetime), fileNumber, length(files), signalNumber, signalNumber + chunkChannelCount-1);
          if(~skipComputation)
              data = load(fullfile(fourierPath, files(fileNumber).name));
              chunkStart = (fileNumber-1) * length(testData.times) + 1;
              chunkEnd = chunkStart + length(data.times) - 1;
          end

          finishedPower = false;
          while(chunkChannelIndex <= chunkChannelCount && ~finishedPower)
              if(~skipComputation)
                  toCopy = data.channels(channelNumber).absolutePower(:, freqNumber);
                  chunk(chunkChannelIndex, chunkStart:chunkEnd) = toCopy';
              end
              chunkChannelIndex = chunkChannelIndex + 1;
              freqNumber = freqNumber + 1;
              if(freqNumber > freqCount)
                  freqNumber = 1;
                  channelNumber = channelNumber + 1;
                  if(channelNumber > length(testData.channels))
                      finishedPower = true;
                  end
              end
          end
          
          %           endIndex = sampleIndex+length(piece) - 1;
          %           fullSignal(sampleIndex:endIndex) = piece;
          %           fullTimes(sampleIndex:endIndex) = data.times;
          if(~skipComputation)
          fullTimes(chunkStart:chunkEnd) = data.times';
          end
          %           sampleIndex = sampleIndex + length(piece);
      end
      if(~skipComputation)
      removeCols = all(isnan(chunk),1);
      chunk(:, removeCols) = [];
      fullTimes(removeCols) = [];
      removeRows = all(isnan(chunk),2);
      chunk(removeRows,:) = [];      
      
      %convolve the loaded chunk
      chunkChannelIndex = oldChunkChannelIndex;
      freqNumber = oldFreqNumber;
      channelNumber = oldChanNumber;      
      end
      while(channelNumber <= length(testData.channels) && chunkChannelIndex <= size(chunk, 1))
          measureLabel = 'absolutePower';
          chanLabel = testData.channels(channelNumber).label;
          %for freqNumber = 1:freqCount
          while(freqNumber <= freqCount && chunkChannelIndex <= size(chunk, 1))
              freqLabel = sprintf('%d-%dHz', testData.freqInfo.lowFrequencies(freqNumber), testData.freqInfo.highFrequencies(freqNumber));
              label = sprintf('%s %s %s', measureLabel, chanLabel, freqLabel);
              fprintf('%s convolving %d of %d (%s)\n', char(datetime), signalNumber, signalCount, label);
              %     sig = fourierEeg.channels(channelNumber).absolutePower(:, freqNumber);
              %     cSig = s7convolveSignal(sig, sampleRate, downsampleRate, markerTimes, fourierEeg.x);
              if(~skipComputation)
                  cSig = s7convolveSignal(chunk(chunkChannelIndex, :), sampleRate, downsampleRate, markerTimes, fullTimes);
                  fprintf('%s done\n', char(datetime));
              end
              
              convolvedFourier.labels{signalNumber} = label;
              convolvedFourier.signals(signalNumber, :) = cSig;
              signalNumber = signalNumber + 1;
              freqNumber = freqNumber + 1;
              chunkChannelIndex = chunkChannelIndex + 1;
          end
          if(chunkChannelIndex <= size(chunk, 1) || freqNumber > freqCount)
              freqNumber = 1;
              channelNumber = channelNumber + 1;
          end
      end      
      oldChunkChannelIndex = 1;
      chunkChannelIndex = 1;
      oldFreqNumber = freqNumber;
      oldChanNumber = channelNumber;
  end
end



%%coherence

chanPairNumber = 1;
chunkChannelIndex = 1;
oldChunkChannelIndex = 1;
while(chanPairNumber <= length(testData.channelPairs))
  measureLabel = 'coherence';
  chanLabel = testData.channelPairs(chanPairNumber).label;
  %for freqNumber = 1:freqCount
  freqNumber = 1;
  while(freqNumber <= freqCount)
      %       fullSignal = NaN(totalSamples, 1);
      fullTimes = NaN(totalSamples, 1);
      %       sampleIndex = 1;
      oldChunkChannelIndex = chunkChannelIndex;
      oldFreqNumber = freqNumber;
      oldChanNumber = chanPairNumber;
      fprintf('%s convolving coherence, chunk %d-%d of %d\n', char(datetime), signalNumber, signalNumber + chunkChannelCount-1, signalCount);
      
      %to conserve ram, copy saved data to array in chunks
      for fileNumber = 1:length(files )

          %piece = data.channelPairs(chanPairNumber).coherence(:, freqNumber);
          chunkChannelIndex = oldChunkChannelIndex;
          freqNumber = oldFreqNumber;
          chanPairNumber = oldChanNumber;
           %           fprintf('%s file %d of %d, chunk %d\n', char(datetime), fileNumber, length(files), chunkChannelIndex);
          fprintf('%s file %d of %d, chunk %d-%d\n', char(datetime), fileNumber, length(files), signalNumber, signalNumber + chunkChannelCount-1);
          data = load(fullfile(fourierPath, files(fileNumber).name));
          chunkStart = (fileNumber-1) * length(testData.times) + 1;
          chunkEnd = chunkStart + length(data.times) - 1;

          finishedPower = false;
          while(chunkChannelIndex <= chunkChannelCount && ~finishedPower)
              toCopy = data.channelPairs(chanPairNumber).coherence(:, freqNumber);
              chunk(chunkChannelIndex, chunkStart:chunkEnd) = toCopy';
              chunkChannelIndex = chunkChannelIndex + 1;
              freqNumber = freqNumber + 1;
              if(freqNumber > freqCount)
                  freqNumber = 1;
                  chanPairNumber = chanPairNumber + 1;
                  if(chanPairNumber > length(testData.channelPairs))
                      finishedPower = true;
                  end
              end
          end
          
          %           endIndex = sampleIndex+length(piece) - 1;
          %           fullSignal(sampleIndex:endIndex) = piece;
          %           fullTimes(sampleIndex:endIndex) = data.times;
          fullTimes(chunkStart:chunkEnd) = data.times';
          %           sampleIndex = sampleIndex + length(piece);
      end
      removeCols = all(isnan(chunk),1);
      chunk(:, removeCols) = [];
      fullTimes(removeCols) = [];
      removeRows = all(isnan(chunk),2);
      chunk(removeRows,:) = [];
      
      %convolve the loaded chunk
      chunkChannelIndex = oldChunkChannelIndex;
      freqNumber = oldFreqNumber;
      chanPairNumber = oldChanNumber;
      while(chanPairNumber <= length(testData.channelPairs) && chunkChannelIndex <= size(chunk, 1))
          measureLabel = 'coherence';
          chanLabel = testData.channelPairs(chanPairNumber).label;
          %for freqNumber = 1:freqCount
          
          while(freqNumber <= freqCount && chunkChannelIndex <= size(chunk, 1))
              freqLabel = sprintf('%d-%dHz', data.freqInfo.lowFrequencies(freqNumber), data.freqInfo.highFrequencies(freqNumber));
              label = sprintf('%s %s %s', measureLabel, chanLabel, freqLabel);
              fprintf('%s convolving %d of %d (%s)\n', char(datetime), signalNumber, signalCount, label);
              %     sig = fourierEeg.channelPairs(chanPairNumber).coherence(:, freqNumber);
              %     cSig = s7convolveSignal(sig, sampleRate, downsampleRate, markerTimes, fourierEeg.x);
              cSig = s7convolveSignal(chunk(chunkChannelIndex, :), sampleRate, downsampleRate, markerTimes, fullTimes);
              fprintf('%s done\n', char(datetime));
              
              convolvedFourier.labels{signalNumber} = label;
              convolvedFourier.signals(signalNumber, :) = cSig;
              signalNumber = signalNumber + 1;
              freqNumber = freqNumber + 1;
              chunkChannelIndex = chunkChannelIndex + 1;
              
          end
          if(chunkChannelIndex <= size(chunk, 1) || freqNumber > freqCount)
              chanPairNumber = chanPairNumber + 1;
              freqNumber = 1;
          end
      end      
      oldChunkChannelIndex = 1;
      chunkChannelIndex = 1;
      oldFreqNumber = freqNumber;
      oldChanNumber = chanPairNumber;
  end
end


timeToFinish = toc;

timeToFinish 
dummy = 1;


% 
% for chanPairNumber = 1:length(data.channelPairs)
%   measureLabel = 'coherence';
%   chanLabel = data.channelPairs(chanPairNumber).label;
%   for freqNumber = 1:freqCount
%             %added
%       fprintf('%s loading %d of %d\n', char(datetime), signalNumber, signalCount);
%       fullSignal = NaN(totalSamples, 1);
%       fullTimes = NaN(totalSamples, 1);
%       sampleIndex = 1;
%       for(fileNumber = 1:length(files ))
%           data = load(fullfile(fourierPath, files(fileNumber).name));
%           piece = data.channelPairs(chanPairNumber).coherence(:, freqNumber);
%           endIndex = sampleIndex+length(piece) - 1;
%           fullSignal(sampleIndex:endIndex) = piece;
%           fullTimes(sampleIndex:endIndex) = data.times;
%           sampleIndex = sampleIndex + length(piece);
%       end
%       removeCols = isnan(fullSignal);
%       fullSignal(removeCols) = [];
%       fullTimes(removeCols) = [];
%       %end added
% 
%     freqLabel = sprintf('%d-%dHz', data.freqInfo.lowFrequencies(freqNumber), data.freqInfo.highFrequencies(freqNumber));
%     label = sprintf('%s %s %s', measureLabel, chanLabel, freqLabel);
%     fprintf('%s convolving %d of %d (%s)\n', char(datetime), signalNumber, signalCount, label);
% %     sig = fourierEeg.coh(chanPairNumber).coherence(:, freqNumber);
% %     cSig = s7convolveSignal(sig, sampleRate, downsampleRate, markerTimes, fourierEeg.x);
%     cSig = s7convolveSignal(fullSignal, sampleRate, downsampleRate, markerTimes, fullTimes);
%     convolvedFourier.labels{signalNumber} = label;
%     convolvedFourier.signals(signalNumber, :) = cSig;
%     signalNumber = signalNumber + 1;
%   end
% end


