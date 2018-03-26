function [ coherencePlot, absolutePowerPlot, phaseAnglePlot ] = coherence( data1, data2, lowFreq, highFreq)
%COHERENCE Computes a timeseries of coherence for a given frequency.
%   Detailed explanation goes here

downSample = false;

sampleRate = 2048;
fftWindowDuration = 64;
refreshRate = 128;

if(length(data1)~=length(data2))
    error('data arrays must be the same length.');
end
if(downSample)
    downSampleFactor = sampleRate / refreshRate;
    newLength= floor(length(data1) / downSampleFactor);
    tempD1 = NaN(1,newLength);
    tempD2 = NaN(1,newLength);
    sourceCounter = 1;
    for i = 1:newLength
        sourceEnd = sourceCounter+downSampleFactor-1;
        tempD1(i) = mean(data1(sourceCounter:sourceEnd));
        tempD2(i) = mean(data2(sourceCounter:sourceEnd));
        sourceCounter = sourceCounter + downSampleFactor;
    end
    sampleRate = refreshRate;
    data1 = tempD1;
    data2 = tempD2;
end


%may want to make these adjustable
if(~exist('lowFreq','var'))
    lowFreq = [1 4 8 12 25 .5 .25 .125];
    highFreq = [4 8 12 25 30 1 .5 .25];
end
coherenceMemoryDuration = 0.5;
if(~exist('lowFreq', 'var'))
    lowFreq = 1;
    processFrequenciesInGroups = false;
else
    processFrequenciesInGroups = true;
end
if(~exist('highFreq', 'var'))
    highFreq = 40;
end
debugMode = false;


%initialize some vars
if(~exist('data1', 'var'))
    if(debugMode)
        simSize = 30 * 2048;
        sampleFreq = 2;
        data1 = sin((1:simSize).*(2*pi*sampleFreq/sampleRate));
        data2 = sin((1:simSize).*(2*pi*sampleFreq/sampleRate) + 2);
        data2(1:(simSize/2)) = sin((1:(simSize/2)).*(2*pi*sampleFreq/sampleRate));
        if(debugMode)
            %debug
            close all;
            figure;
            hold on;
            plot(data1);
            plot(data2,'r');
            %end debug
        end
    else
        error('No data was passed to this function');
    end
end
plotLength=floor(length(data1) / (sampleRate/refreshRate) - refreshRate * fftWindowDuration);
if(~processFrequenciesInGroups)
    freqCount = highFreq-lowFreq+1;
else
    freqCount = length(lowFreq);
end
coherencePlot = NaN(plotLength, freqCount);
absolutePowerPlot = NaN(plotLength, freqCount);
if(nargout > 2)
    phaseAnglePlot = NaN(plotLength, freqCount);
end

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
fftFrequencies = (1:fftEndIndex) ./ fftWindowDuration;

% tickSize = ceiling(length(data1 / 256));
if(~processFrequenciesInGroups)
    iRange = (lowFreq+1):(highFreq+1);
else
    iRange = 1:length(lowFreq);
end

while(fftEndIndex < length(data1))
    fft1 = fft(data1(fftStartIndex:fftEndIndex));
    fft2 = fft(data2(fftStartIndex:fftEndIndex));
    if(processFrequenciesInGroups)
        tempFft1 = complex(NaN(1,length(iRange)));
        tempFft2 = complex(NaN(1,length(iRange)));
        for i = 1:length(lowFreq)
            tempFft1(i) = mean(fft1((lowFreq(i)*fftWindowDuration+1):highFreq(i)*fftWindowDuration+1));
            tempFft2(i) = mean(fft2((lowFreq(i)*fftWindowDuration+1):highFreq(i)*fftWindowDuration+1));
        end
        fft1 = tempFft1;
        fft2 = tempFft2;
    end
    for i = iRange
        if(coherenceFull)
            powerSum1(i) = powerSum1(i) - powerArray1(coherenceIndex, i);
            powerSum2(i) = powerSum2(i) - powerArray2(coherenceIndex, i);
            dotProdSum(i) = dotProdSum(i) - dotProdArray(coherenceIndex, i);
            crossProdSum(i) = crossProdSum(i) - crossProdArray(coherenceIndex, i);
        end
        x1 = real(fft1(i));
        y1 = imag(fft1(i));
        x2 = real(fft2(i));
        y2 = imag(fft2(i));
        %multiply by 2 because you'll be adding dot and cross products together
        pow1 = 2*(x1 * x1 + y1 * y1);
        pow2 = 2*(x2 * x2 + y2 * y2);
        dot = 2*(x1 * x2 + y1 * y2);
        cross = 2*(x1 * y2 - x2 * y1);
        
        powerArray1(coherenceIndex, i) = pow1;
        powerArray2(coherenceIndex, i) = pow2;
        dotProdArray(coherenceIndex, i) = dot;
        crossProdArray(coherenceIndex, i) = cross;
        
        powerSum1(i) = powerSum1(i) + powerArray1(coherenceIndex, i);
        powerSum2(i) = powerSum2(i) + powerArray2(coherenceIndex, i);
        dotProdSum(i) = dotProdSum(i) + dotProdArray(coherenceIndex, i);
        crossProdSum(i) = crossProdSum(i) + crossProdArray(coherenceIndex, i);
        if(coherenceFull)
            %todo: do we need to divide by queueLength here?
            avgPow1 = powerSum1(i) / queueLength;
            avgPow2 = powerSum2(i) / queueLength;
            avgDot = dotProdSum(i) / queueLength;
            avgCross = crossProdSum(i) / queueLength;
            coherencePlot(coherenceCount, i) = (avgDot*avgDot + avgCross*avgCross)...
                / (avgPow1*avgPow2);
            absolutePowerPlot(coherenceCount, i) = (avgPow1);
            if(nargout > 1)
                for phaseCounter = 1:freqCount
                    myPhase = phaseAngle(fft1,fft2,lowFreq,highFreq);
                    phaseAnglePlot(coherenceCount, phaseCounter) = myPhase;
                end
            end
        end
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

if(coherenceCount <=size(coherencePlot, 1))
    phaseAnglePlot(coherenceCount:size(coherencePlot,1),:) = [];
    coherencePlot(coherenceCount:size(coherencePlot,1),:) = [];
end

if(debugMode)
    %debug
    for i = iRange
        close all;
        plot(coherencePlot(:,i));
        if(~processFrequenciesInGroups)
            title(sprintf('%d', i-1));
        else
            title(sprintf('%d - %d Hz', lowFreq(i), highFreq(i)));
        end
    end
    %end debug
end

end

