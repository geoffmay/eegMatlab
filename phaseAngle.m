function [ phaseDiff, x1, y1, x2, y2 ] = phaseAngle( fft1, fft2, lowerFreqIndex, upperFreqIndex )
%PHASEANGLE Computes the phase angle between two fast fourier transforms at
%a given range of frequency indexes.
%   Returns phase diff, and optionally, the summed components for power
%   calculation.

if(lowerFreqIndex == 1)
    warning('The first frequency index is the dc offset.  Did you mean to start with 2?');
end

%grab the frequency range from the fft
ff1 = fft1(lowerFreqIndex:upperFreqIndex);
ff2 = fft2(lowerFreqIndex:upperFreqIndex);

%add the vectors, weighting each frequency equally.
x1 = sum(real(ff1));
y1 = sum(imag(ff1));
x2 = sum(real(ff2));
y2 = sum(imag(ff2));

%compute phase angles
phase1 = atan2(y1, x1);
phase2 = atan2(y2, x2);

phaseDiff = phase2 - phase1;

if(phaseDiff < -pi)
    phaseDiff = 2*pi + phaseDiff;
elseif(phaseDiff > pi)
    phaseDiff = -2*pi + phaseDiff;
end


end

