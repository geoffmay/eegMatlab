function [ slope, meanPhaseAngle, polyLine ] = phaseSlopeIndexOld( data, chanIndex1, chanIndex2 )
%PHASESLOPEINDEX Summary of this function goes here
%   Detailed explanation goes here

sampleRate = 2048;
windowIncrement = 128;
subSampleStartFreq = 20;
subSampleEndFreq = 100;
sampleCounter = 1;
totalWindows = floor(size(data,1) / sampleRate);
windowCounter = 1;
phaseAngles = NaN(totalWindows, windowIncrement);

while(sampleCounter + sampleRate < size(data,1))
    endIndex = sampleCounter + sampleRate - 1;
    sample1 = data(sampleCounter:endIndex, chanIndex1);
    sample2 = data(sampleCounter:endIndex, chanIndex2);
    fft1 = fft(sample1);
    fft2 = fft(sample2);
    for i = 1:windowIncrement
        phase1 = atan2(imag(fft1(i)), real(fft1(i)));
        phase2 = atan2(imag(fft2(i)), real(fft2(i)));
        if(phase1 < pi)
            phase1 = phase1 + 2*pi;
        elseif(phase1 > pi)
            phase1 = phase1 - 2*pi;
        end
        if(phase2 < pi)
            phase2 = phase2 + 2*pi;
        elseif(phase2 > pi)
            phase2 = phase2 - 2*pi;
        end
        phaseAngles(windowCounter, i) = phase2 - phase1;
    end
    sampleCounter = sampleCounter + windowIncrement;
    windowCounter = windowCounter + 1;
end
meanPhaseAngle = mean(phaseAngles,1);
[polyLine, error] = polyfit(subSampleStartFreq:subSampleEndFreq, meanPhaseAngle(subSampleStartFreq:subSampleEndFreq), 1);
slope = polyLine(1);



end

