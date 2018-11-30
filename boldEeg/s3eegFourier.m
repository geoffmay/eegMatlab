function s3eegFourier( EEG, outputFolder )
%COHERENCE Computes a timeseries of coherence for a given frequency.
%   Detailed explanation goes here


highRes = 0;

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
    lowFreq =  [1 5 9  13 25 40 70];
    highFreq = [4 8 12 24 39 70 124];
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

if(EEG.srate == 250)
    refreshRate = 250;
elseif(EEG.srate == 1000)
    refreshRate = 250;
end

freqInfo.lowFrequencies = lowFreq;
freqInfo.highFrequencies = highFreq;
freqInfo.coherenceSampleRateHz = refreshRate;
freqInfo.coherenceMemoryDurationSeconds = coherenceMemoryDuration;
freqInfo.fftWindowDurationSeconds = fftWindowDuration;
freqInfo.samplesPerChunk = refreshRate * 10;
freqInfo.chunkCounter = 1;
freqInfo.sourceFile = EEG.filepath;



if(~exist('EEG', 'var'))
  error('No data was passed to this function');
end
mkdir(outputFolder);

%initialize some vars
% plotLength=floor(size(EEG.data,2) / (EEG.srate/refreshRate) - refreshRate * fftWindowDuration);
plotLength = freqInfo.samplesPerChunk;
channelPairCount = EEG.nbchan * (EEG.nbchan-1) / 2;
nullPlot = NaN(plotLength,length(freqInfo.lowFrequencies));
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


%initialize arrays
queueLength = refreshRate * coherenceMemoryDuration;
queueLengthRecip = 1/queueLength;
powerArray = NaN(queueLength, freqCount, EEG.nbchan);
powerSum = zeros(freqCount, EEG.nbchan);
powerAvg = zeros(freqCount, EEG.nbchan);
dotProdArray = NaN(queueLength, freqCount, channelPairCount);
crossProdArray = NaN(queueLength, freqCount, channelPairCount);
dotProdSum = zeros(freqCount, channelPairCount);
crossProdSum = zeros(freqCount, channelPairCount);

coherenceIndex = 1;
coherenceCount = 1;
coherenceFull = false;
chunkIndex = 1;
fftStartIndex = 1;
fftEndIndex = fftStartIndex + sampleRate*fftWindowDuration - 1;
ffts = NaN(fftEndIndex, channelPairCount);
meanFfts = NaN(length(lowFreq), EEG.nbchan);
freqRange = 1:length(lowFreq);

fprintf('-');


while(fftEndIndex < size(EEG.data,2))
  times(chunkIndex) = coherenceCount / refreshRate;
  fprintf('\b');
  ticker = mod(coherenceCount, 4);
  if(ticker==0)
    fprintf('\\');
  elseif(ticker==1)
    fprintf('-');
  elseif(ticker==2)
    fprintf('/');
  elseif(ticker==3)
    fprintf('|');
  end
  if(mod(coherenceCount, 1000)==0)
    percentDone =  100 * fftEndIndex / size(EEG.data,2);
    fprintf('\n%s (%.3f%%)-', char(datetime), percentDone);
  end
  for chanCounter = 1:EEG.nbchan
%     tic;
    ffts(:, chanCounter) = fft(EEG.data(chanCounter, fftStartIndex:fftEndIndex));
%     clock.fft = toc;
%     tic;
    for i = iRange
      if(processFrequenciesInGroups)
        myRange = (lowFreq(i)*fftWindowDuration+1):(highFreq(i)*fftWindowDuration+1);
        meanFfts(i, chanCounter) = mean(ffts(myRange, chanCounter));
      else
        meanFfts(i, chanCounter) = ffts(i, chanCounter);
        
      end
      if(coherenceFull)
        powerSum(i, chanCounter) = powerSum(i, chanCounter) - powerArray(coherenceIndex, i, chanCounter);
      end
      x = real(meanFfts(i,chanCounter));
      y = imag(meanFfts(i,chanCounter));
      pow = 2*(x * x + y * y);
      powerArray(coherenceIndex, i, chanCounter) = pow;
      powerSum(i, chanCounter) = powerSum(i, chanCounter) + powerArray(coherenceIndex, i, chanCounter);
      avgPow(i,chanCounter) = powerSum(i, chanCounter) * queueLengthRecip;
    end
%     clock.average = toc;

  end

  
  channelPairIndex = 1;
  for chan1 = 1:EEG.nbchan
    for chan2 = chan1+1:EEG.nbchan
      if(length(chan2) > 0)
        for i = iRange
%           tic;
          if(coherenceFull)
            dotProdSum(i, channelPairIndex) = dotProdSum(i, channelPairIndex) ...
              - dotProdArray(coherenceIndex, i, channelPairIndex);
            crossProdSum(i, channelPairIndex) = crossProdSum(i, channelPairIndex) ...
              - crossProdArray(coherenceIndex, i, channelPairIndex);
          end
          x1 = real(meanFfts(i,chan1));
          y1 = imag(meanFfts(i,chan1));
          x2 = real(meanFfts(i,chan2));
          y2 = imag(meanFfts(i,chan2));
          dot = 2*(x1 * x2 + y1 * y2);
          cross = 2*(x1 * y2 - x2 * y1);
          dotProdArray(coherenceIndex, i, channelPairIndex) = dot;
          crossProdArray(coherenceIndex, i, channelPairIndex) = cross;
          
          dotProdSum(i, channelPairIndex) = dotProdSum(i, channelPairIndex)...
            + dotProdArray(coherenceIndex, i, channelPairIndex);
          crossProdSum(i, channelPairIndex) = crossProdSum(i, channelPairIndex)...
            + crossProdArray(coherenceIndex, i, channelPairIndex);
          if(coherenceFull)
            avgPow1 = avgPow(i, chan1);
            avgPow2 = avgPow(i, chan2);
            avgDot = dotProdSum(i, channelPairIndex) * queueLengthRecip;
            avgCross = crossProdSum(i, channelPairIndex) * queueLengthRecip;
            channelPairs(channelPairIndex).coherence(chunkIndex, i) = ...
              (avgDot*avgDot + avgCross*avgCross)...
              / (avgPow1*avgPow2);
            channels(chan1).absolutePower(chunkIndex, i) = log(avgPow1);
            %need to fill in the last channel...
            if(chan1 == 1)
              channels(chan2).absolutePower(chunkIndex, i) = log(avgPow2);
            end
%           clock.coherence = toc;

          end
        end
        channelPairIndex = channelPairIndex + 1;
      end
    end
    if(coherenceFull)
      thisPowerSumRecip = 1 / sum(channels(chan1).absolutePower(chunkIndex, :));
      channels(chan1).relativePower(chunkIndex, :) = channels(chan1).absolutePower(chunkIndex, :) .* thisPowerSumRecip;
    end
  end
  %     %debug
  %     if(coherenceFull)
  %         avgPow1 = avgPow(i, chan1);
  %     end
  %     %end debug
  if(chunkIndex == freqInfo.samplesPerChunk)
      %write file
      outputFilename = fullfile(outputFolder, sprintf('%05d.mat', freqInfo.chunkCounter));
      save(outputFilename, 'times', 'channels', 'channelPairs', 'freqInfo', '-v7.3');
      freqInfo.chunkCounter = freqInfo.chunkCounter + 1;
      chunkIndex = 1;
  end

  chunkIndex = mod(coherenceCount, freqInfo.samplesPerChunk) + 1;
  coherenceCount = coherenceCount + 1;
  
  coherenceIndex = coherenceIndex + 1;
  if(coherenceIndex > queueLength)
    coherenceIndex = 1;
    coherenceFull = true;
  end
  fftStartIndex = fftStartIndex + (sampleRate / refreshRate);
  fftEndIndex = fftStartIndex + sampleRate*fftWindowDuration - 1;
end

if(chunkIndex > 0)
    if(chunkIndex < length(times))
        times((chunkIndex + 1):end) = [];
        for i = 1:length(channelPairs)
            channelPairs(i).coherence((chunkIndex + 1):end,:) = [];
        end
        for i = 1:length(channels)
            channels(i).absolutePower((chunkIndex + 1):end,:)=[];
            channels(i).relativePower((chunkIndex + 1):end,:)=[];
        end
        
    end
    outputFilename = fullfile(outputFolder, sprintf('%05d.mat', freqInfo.chunkCounter));
    save(outputFilename, 'times', 'channels', 'channelPairs', 'freqInfo');
end

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

