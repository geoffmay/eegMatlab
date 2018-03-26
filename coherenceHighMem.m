function [ coherencePlot, absolutePowerPlot, phaseAnglePlot ] = coherence( allChannelData, lowFreq, highFreq)
%COHERENCE Computes a timeseries of coherence for a given frequency.
%   Detailed explanation goes here

%flags
downSample = false;
debugMode = false;

%adjustable constants
sampleRate = 2048;
fftWindowDuration = 64;
refreshRate = 128;
coherenceMemoryDuration = 0.5;


if(~exist('lowFreq','var'))
    lowFreq = [1 4 8 12 25 .5 .25 .125];
    highFreq = [4 8 12 25 30 1 .5 .25];
    %     lowFreq = [1 4 8 12];
    %     highFreq = [4 8 12 25];
end

%derive sizes
freqCount = length(lowFreq);
plotLength=floor(size(allChannelData, 1) / (sampleRate/refreshRate) - refreshRate * fftWindowDuration);
channelCount = size(allChannelData, 2);
% %debug
% plotLength = 1000;
% channelCount = 4;
% %end debug

queueLength = refreshRate * coherenceMemoryDuration;
maxFrequency= max(highFreq);
fourierSize = maxFrequency * fftWindowDuration;

%initialize arrays
coherencePlot = NaN(plotLength, freqCount, channelCount);
absolutePowerPlot = NaN(plotLength, freqCount, channelCount);
phaseAnglePlot = NaN(plotLength, freqCount, channelCount);

ffts = complex(NaN(plotLength, freqCount, channelCount));
% powerArray1 = NaN(queueLength, freqCount);
% powerArray2 = NaN(queueLength, freqCount);
% powerSum1 = zeros(freqCount, 1);
% powerSum2 = zeros(freqCount, 1);
% dotProdArray = NaN(queueLength, freqCount);
% crossProdArray = NaN(queueLength, freqCount);
% dotProdSum = zeros(freqCount, 1);
% crossProdSum = zeros(freqCount, 1);

fftStartIndex = 1;
fftEndIndex = fftStartIndex + sampleRate*fftWindowDuration - 1;
fftFrequencies = (1:fftEndIndex) ./ fftWindowDuration;

%compute fourier transforms
for chanCounter = 1:channelCount
    fprintf('\n%s computing transform %d of %d', char(datetime), chanCounter, channelCount);
    for freqCounter = 1:freqCount
        fprintf('.');
        fftStartIndex = 1;
        fftEndIndex = fftStartIndex + sampleRate*fftWindowDuration - 1;
        plotCounter = 1;
        while(plotCounter < plotLength)
            fft1 = fft(allChannelData(fftStartIndex:fftEndIndex, chanCounter));
            fft1 = log(fft1);
            tempFft1 = complex(NaN(1,freqCount));
            for i = 1:freqCount
                tempFft1(i) = mean(fft1((lowFreq(i)*fftWindowDuration+1):highFreq(i)*fftWindowDuration+1));
            end
            fft1 = tempFft1;
            ffts(plotCounter,:,chanCounter) = fft1;
            fftStartIndex = fftStartIndex + (sampleRate / refreshRate);
            fftEndIndex = fftStartIndex + sampleRate*fftWindowDuration - 1;
            plotCounter = plotCounter + 1;
        end
    end
end

%figure out how many channel pairs there are
maxPairs = 0;
indexTable = [];
%pairLabels = cell(0);
for chan1 = 1:(channelCount)
    for chan2 = chan1+1:(channelCount)
        maxPairs = maxPairs + 1;
        indexTable(maxPairs, 1) = maxPairs;
        indexTable(maxPairs, 2) = chan1;
        indexTable(maxPairs, 3) = chan2;
        %pairLabels{maxPairs}= sprintf('%s-%s',labels{chan1},labels{chan2});
    end
end

%compute coherence
chanPairCounter = 1;
for chan1 = 1:channelCount
    for chan2 = chan1+1:channelCount
        for freqCounter = 1:freqCount
            fprintf('\ncoherence %d of %d', chanPairCounter, maxPairs);            
            coherenceIndex = 1;
            coherenceCounter = 1;
            coherenceFull = false;
            powerArray1 = NaN(queueLength, 1);
            powerArray2 = NaN(queueLength, 1);
            dotProdArray = NaN(queueLength, 1);
            crossProdArray = NaN(queueLength, 1);
            powerSum1 = 0;
            powerSum2 = 0;
            dotProdSum = 0;
            crossProdSum = 0;
            for plotCounter = 1:size(ffts,1)                
                if(coherenceFull)
                    powerSum1 = powerSum1 - powerArray1(coherenceIndex);
                    powerSum2 = powerSum2 - powerArray2(coherenceIndex);
                    dotProdSum = dotProdSum - dotProdArray(coherenceIndex);
                    crossProdSum = crossProdSum - crossProdArray(coherenceIndex);
                end
                x1 = real(ffts(plotCounter, freqCounter, chan1));
                y1 = imag(ffts(plotCounter, freqCounter, chan1));
                x2 = real(ffts(plotCounter, freqCounter, chan2));
                y2 = imag(ffts(plotCounter, freqCounter, chan2));
                %multiply by 2 because you'll be adding dot and cross products together
                pow1 = 2*(x1 * x1 + y1 * y1);
                pow2 = 2*(x2 * x2 + y2 * y2);
                dot = 2*(x1 * x2 + y1 * y2);
                cross = 2*(x1 * y2 - x2 * y1);
                
                powerArray1(coherenceIndex) = pow1;
                powerArray2(coherenceIndex) = pow2;
                dotProdArray(coherenceIndex) = dot;
                crossProdArray(coherenceIndex) = cross;
                
                powerSum1 = powerSum1 + powerArray1(coherenceIndex);
                powerSum2 = powerSum2 + powerArray2(coherenceIndex);
                dotProdSum = dotProdSum + dotProdArray(coherenceIndex);
                crossProdSum = crossProdSum + crossProdArray(coherenceIndex);
                if(coherenceFull)
                    %todo: do we need to divide by queueLength here?
                    avgPow1 = powerSum1 / queueLength;
                    avgPow2 = powerSum2 / queueLength;
                    avgDot = dotProdSum / queueLength;
                    avgCross = crossProdSum / queueLength;
                    coherencePlot(coherenceCounter, freqCounter, chanPairCounter) = (avgDot*avgDot + avgCross*avgCross)...
                        / (avgPow1*avgPow2);
                    absolutePowerPlot(coherenceCounter, i) = (avgPow1);
                    if(false)
                        for phaseCounter = 1:freqCount
                            myPhase = phaseAngle(fft1,fft2,lowFreq,highFreq);
                            phaseAnglePlot(coherenceCounter, phaseCounter) = myPhase;
                        end
                    end
                end
                coherenceIndex = coherenceIndex + 1;
                if(coherenceFull)
                    coherenceCounter = coherenceCounter + 1;
                end
                if(coherenceIndex > queueLength)
                    coherenceIndex = 1;
                    coherenceFull = true;
                end
            end
        end
        chanPairCounter = chanPairCounter + 1;
    end
end

% 
% %old code
% 
% iRange = 1:length(lowFreq);
% while(fftEndIndex < length(data1))
%     fft1 = fft(data1(fftStartIndex:fftEndIndex));
%     fft2 = fft(data2(fftStartIndex:fftEndIndex));
%     tempFft1 = complex(NaN(1,length(iRange)));
%     tempFft2 = complex(NaN(1,length(iRange)));
%     for i = 1:length(lowFreq)
%         tempFft1(i) = mean(fft1((lowFreq(i)*fftWindowDuration+1):highFreq(i)*fftWindowDuration+1));
%         tempFft2(i) = mean(fft2((lowFreq(i)*fftWindowDuration+1):highFreq(i)*fftWindowDuration+1));
%     end
%     fft1 = tempFft1;
%     fft2 = tempFft2;
%     for i = iRange
%         if(coherenceFull)
%             powerSum1(i) = powerSum1(i) - powerArray1(coherenceIndex, i);
%             powerSum2(i) = powerSum2(i) - powerArray2(coherenceIndex, i);
%             dotProdSum(i) = dotProdSum(i) - dotProdArray(coherenceIndex, i);
%             crossProdSum(i) = crossProdSum(i) - crossProdArray(coherenceIndex, i);
%         end
%         x1 = real(fft1(i));
%         y1 = imag(fft1(i));
%         x2 = real(fft2(i));
%         y2 = imag(fft2(i));
%         %multiply by 2 because you'll be adding dot and cross products together
%         pow1 = 2*(x1 * x1 + y1 * y1);
%         pow2 = 2*(x2 * x2 + y2 * y2);
%         dot = 2*(x1 * x2 + y1 * y2);
%         cross = 2*(x1 * y2 - x2 * y1);
%         
%         powerArray1(coherenceIndex, i) = pow1;
%         powerArray2(coherenceIndex, i) = pow2;
%         dotProdArray(coherenceIndex, i) = dot;
%         crossProdArray(coherenceIndex, i) = cross;
%         
%         powerSum1(i) = powerSum1(i) + powerArray1(coherenceIndex, i);
%         powerSum2(i) = powerSum2(i) + powerArray2(coherenceIndex, i);
%         dotProdSum(i) = dotProdSum(i) + dotProdArray(coherenceIndex, i);
%         crossProdSum(i) = crossProdSum(i) + crossProdArray(coherenceIndex, i);
%         if(coherenceFull)
%             %todo: do we need to divide by queueLength here?
%             avgPow1 = powerSum1(i) / queueLength;
%             avgPow2 = powerSum2(i) / queueLength;
%             avgDot = dotProdSum(i) / queueLength;
%             avgCross = crossProdSum(i) / queueLength;
%             coherencePlot(coherenceCounter, i) = (avgDot*avgDot + avgCross*avgCross)...
%                 / (avgPow1*avgPow2);
%             absolutePowerPlot(coherenceCounter, i) = (avgPow1);
%             if(nargout > 1)
%                 for phaseCounter = 1:freqCount
%                     myPhase = phaseAngle(fft1,fft2,lowFreq,highFreq);
%                     phaseAnglePlot(coherenceCounter, phaseCounter) = myPhase;
%                 end
%             end
%         end
%     end
%     
%     coherenceIndex = coherenceIndex + 1;
%     if(coherenceFull)
%         coherenceCounter = coherenceCounter + 1;
%     end
%     if(coherenceIndex > queueLength)
%         coherenceIndex = 1;
%         coherenceFull = true;
%     end
%     fftStartIndex = fftStartIndex + (sampleRate / refreshRate);
%     fftEndIndex = fftStartIndex + sampleRate*fftWindowDuration - 1;
% end

if(coherenceCounter <=size(coherencePlot, 1))
    phaseAnglePlot(coherenceCounter:size(coherencePlot,1),:) = [];
    coherencePlot(coherenceCounter:size(coherencePlot,1),:) = [];
end

% if(debugMode)
%     %debug
%     for i = iRange
%         close all;
%         plot(coherencePlot(:,i));
%         title(sprintf('%d - %d Hz', lowFreq(i), highFreq(i)));
%     end
%     %end debug
% end

end

