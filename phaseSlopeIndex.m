function [ timeLagSeconds, meanPhaseAngle, polyLine, maxTimeShift ] = phaseSlopeIndex( EEG, chanIndex1, chanIndex2, minFrequency, maxFrequency)
%PHASESLOPEINDEX calculates the phase slope between two columns of a matrix
%of timeseries data.
%   Phase slope is an approximation of the way that the phase angle between
%   two timeseries data plots behaves as frequency increases.  If the phase
%   slope is positive, then channel 2 precedes channel 1.  If the phase
%   slope is negative, then channel 1 precedes channel 2 in time.

debug = false;
plot1 = false;
sizeToIncrementRatio = 16;


sampleRate = EEG.srate;
if(floor(minFrequency) ~= minFrequency || floor(maxFrequency) ~= maxFrequency)
    error('currently only integer frequency limits are supported');
end
windowSize = sampleRate;
if(size(EEG.data, 2) < windowSize)
    error('not enough data; need at least %d samples', windowSize);
end
maxFrameShift = floor(EEG.srate/(2*maxFrequency));
%maxTimeShift = 1 / (2*maxFrequency);
maxTimeShift = maxFrameShift / EEG.srate;
windowIncrement = windowSize / sizeToIncrementRatio;
frequenciesOfInterest = minFrequency:maxFrequency;
fftIndices = frequenciesOfInterest + 1;
maxFft = max(fftIndices);
% subSampleStartFreq = 2;
% subSampleEndFreq = 31;
sampleCounter = 1;
totalWindows = floor(size(EEG.data,2) / windowIncrement) - sizeToIncrementRatio;
windowCounter = 1;
phaseAngles = NaN(totalWindows, maxFft);
if(debug)
    allPhase1 = NaN(totalWindows, maxFft);
    allPhase2 = NaN(totalWindows, maxFft);
end
while(sampleCounter + sampleRate < size(EEG.data,2))
    endIndex = sampleCounter + sampleRate - 1;
    sample1 = EEG.data(chanIndex1, sampleCounter:endIndex);
    sample2 = EEG.data(chanIndex2, sampleCounter:endIndex);
    fft1 = fft(sample1);
    fft2 = fft(sample2);
    for freqCounter = fftIndices
        phase1 = atan2(imag(fft1(freqCounter)), real(fft1(freqCounter)));
        phase2 = atan2(imag(fft2(freqCounter)), real(fft2(freqCounter)));
        %positive phase angle means phase2 came first (see debug block for proof)
        phaseAngle = phase2-phase1;
        if(debug)
            close all;
            allPhase1(windowCounter,freqCounter) = phase1;
            allPhase2(windowCounter,freqCounter) = phase2;
            if(plot1)
                if(freqCounter == 3)
                    keepI = [freqCounter+1, length(fft1)-freqCounter+1];
                    keep = zeros(1,length(fft1));
                    keep(keepI) = 1;
                    filtered1 = fft1;
                    filtered2 = fft2;
                    filtered1(find(~keep)) = 0;
                    filtered2(find(~keep)) = 0;
                    fData1 = ifft(filtered1);
                    fData2 = ifft(filtered2);
                    figure;
                    hold on;
                    plot(fData1,'b');
                    plot(fData2,'r');
                end %break here to compare plot to variables.
            end
        end

        %todo: took this out
        %         if(phaseAngle < -pi)
        %             phaseAngle = phaseAngle + 2*pi;
        %         elseif(phaseAngle > pi)
        %             phaseAngle = phaseAngle - 2*pi;
        %         end
        phaseAngles(windowCounter, freqCounter) = phaseAngle;
    end    
    sampleCounter = sampleCounter + windowIncrement;
    windowCounter = windowCounter + 1;
end
if(windowCounter <= totalWindows)
    phaseAngles(windowCounter:end,:) = [];
    if(debug)
        allPhase1(windowCounter:end,:)=[];
        allPhase2(windowCounter:end,:)=[];
    end
end
meanPhaseAngle = mean(phaseAngles,1);

%This step converts phase slope to an estimated time lead (only if window size is 1 second).
meanPhaseAngle = meanPhaseAngle ./ (2*pi);
maxSlopeIndex = min(max(fftIndices),maxFrameShift);

[polyLine, myError] = polyfit(fftIndices, meanPhaseAngle(fftIndices), 1);

%positive slope means that phase 2's slope is greater than phase 1's slope.
%In other words, timecourse 1 lags behind timecourse 2.
timeLagSeconds = polyLine(1);

if(debug)
    meanPhase1 = mean(allPhase1,1);
    meanPhase2 = mean(allPhase2,1);
    close all;
    x = fftIndices;
    [angleFit, myError] = polyfit(x, meanPhaseAngle(x), 1);
    [phase1Fit, myError] = polyfit(x, meanPhase1(x), 1);
    [phase2Fit, myError] = polyfit(x, meanPhase2(x), 1);
    angleY = x .* angleFit(1) + angleFit(2);
    phase1Y = x .* phase1Fit(1) + angleFit(2);
    phase2Y = x .* phase2Fit(1) + angleFit(2);
    figure;
    hold all;
    plot(x, meanPhase1(x), 'b');
    plot(x, meanPhase2(x), 'r');
    plot(x, meanPhaseAngle(x), 'g');
    plot(x, angleY, 'b', 'linewidth', 4);
    plot(x, phase1Y, 'r', 'linewidth', 4);
    plot(x, phase2Y, 'g', 'linewidth', 4);
    legend('phase1', 'phase2', 'phaseAngle', 'p1Fit', 'p2Fit', 'angleFit');
end


end

