function [absPowerPlot, relPowerPlot] = smearedPower( eegData, lowFrequencies, highFrequencies)

downSample = false;
sampleRate = 2048;
fftWindowDuration = 1;
refreshRate = 128;

if(downSample)
  downSampleFactor = sampleRate / refreshRate;
  newLength= floor(length(eegData) / downSampleFactor);
  tempD1 = NaN(1,newLength);
  tempD2 = NaN(1,newLength);
  sourceCounter = 1;
  for i = 1:newLength
    sourceEnd = sourceCounter+downSampleFactor-1;
    tempD1(i) = mean(eegData(sourceCounter:sourceEnd));
    sourceCounter = sourceCounter + downSampleFactor;    
  end
  sampleRate = refreshRate;
  eegData = tempD1;
end


%may want to make these adjustable
if(~exist('lowFreq','var'))
  lowFrequencies = [1 4 8 12 25 30];
  highFrequencies = [4 8 12 25 30 60];
end
coherenceMemoryDuration = 0.5;
if(~exist('lowFrequencies', 'var'))
  lowFrequencies = 1;
  processFrequenciesInGroups = false;
else
  processFrequenciesInGroups = true;
end
if(~exist('highFrequencies', 'var'))
  highFrequencies = 40;
end
debugMode = false;


%initialize some vars
if(~exist('eegData', 'var'))
  if(debugMode)
    simSize = 30 * 2048;
    sampleFreq = 2;
    eegData = sin((1:simSize).*(2*pi*sampleFreq/sampleRate));
    if(debugMode)
      %debug
      close all;
      figure;
      hold on;
      plot(eegData);
      %end debug
    end
  else
    error('No data was passed to this function');
  end
end
plotLength=floor(length(eegData) / (sampleRate/refreshRate) - refreshRate * fftWindowDuration);
if(~processFrequenciesInGroups)
  freqCount = highFrequencies-lowFrequencies+1;
else
  freqCount = length(lowFrequencies);
end
absPowerPlot = NaN(plotLength, freqCount);
relPowerPlot = NaN(plotLength, freqCount);

%initialize arrays
queueLength = refreshRate * coherenceMemoryDuration;
powerArray1 = NaN(queueLength, freqCount);
powerArray2 = NaN(queueLength, freqCount);
powerSum1 = zeros(freqCount, 1);
powerSum2 = zeros(freqCount, 1);
dotProdArray = NaN(queueLength, freqCount);
crossProdArray = NaN(queueLength, freqCount);
dotProdSum = zeros(freqCount, 1);
crossProdSum = zeros(freqCount, 1);

coherenceIndex = 1;
coherenceCount = 1;
coherenceFull = false;
fftStartIndex = 1;
fftEndIndex = fftStartIndex + sampleRate*fftWindowDuration - 1;

% tickSize = ceiling(length(data1 / 256));
if(~processFrequenciesInGroups)
  iRange = (lowFrequencies+1):(highFrequencies+1);
else
  iRange = 1:length(lowFrequencies);
end

while(fftEndIndex < length(eegData))
  fft1 = fft(eegData(fftStartIndex:fftEndIndex));
  if(processFrequenciesInGroups)
    tempFft1 = complex(NaN(1,length(iRange)));
    for i = 1:length(lowFrequencies)
      tempFft1(i) = mean(fft1((lowFrequencies(i)+1):highFrequencies(i)+1));
    end
    fft1 = tempFft1;
  end
  for i = iRange
    if(coherenceFull)
      powerSum1(i) = powerSum1(i) - powerArray1(coherenceIndex, i);
    end
    x1 = real(fft1(i));
    y1 = imag(fft1(i));
    %multiply by 2 because you'll be adding dot and cross products together?
    pow1 = 2*(x1 * x1 + y1 * y1);
    
    powerArray1(coherenceIndex, i) = pow1;
    
    powerSum1(i) = powerSum1(i) + powerArray1(coherenceIndex, i);
    if(coherenceFull)
      %todo: do we need to divide by queueLength here?
      avgPow1 = powerSum1(i) / queueLength;
      absPowerPlot(coherenceCount, i) = log10(avgPow1);

      %       coherencePlot(coherenceCount, i) = (avgDot*avgDot + avgCross*avgCross)...
      %         / (avgPow1*avgPow2);
    end
  end
  sumLogPower = sum(absPowerPlot(coherenceCount, :));
  for i = iRange
      %debug
      if(~any(isnan(sumLogPower)))
          dummy = 1;
      end
      %end debug
      relPowerPlot(coherenceCount, i) = absPowerPlot(coherenceCount,i) / sumLogPower;
  end
  
  coherenceIndex = coherenceIndex + 1;
  if(coherenceFull)
    coherenceCount = coherenceCount + 1;
  end
  if(coherenceIndex > queueLength)
    coherenceIndex = 1;
    coherenceFull = true;
  end
  fftStartIndex = fftStartIndex + (sampleRate / refreshRate);
  fftEndIndex = fftStartIndex + sampleRate*fftWindowDuration - 1;
end

if(coherenceCount <=size(absPowerPlot, 1))
  absPowerPlot(coherenceCount:size(absPowerPlot,1),:) = [];
  relPowerPlot(coherenceCount:size(relPowerPlot,1),:) = [];
end

if(debugMode)
  %debug
  for i = iRange
    close all;
    plot(absPowerPlot(:,i));
    if(~processFrequenciesInGroups)
      title(sprintf('%d', i-1));
    else
      title(sprintf('%d - %d Hz', lowFrequencies(i), highFrequencies(i)));
    end
  end
  %end debug
end


end