function [ filtered ] = viewFilter( gmd, plotTime, windowSize, color )
%VIEWFILTER Summary of this function goes here
%   Detailed explanation goes here

b = (1/windowSize)*ones(1,windowSize);
filt = filtfilt(b,1,gmd);
plot(plotTime, filt, color);

[fPeaks, fPeakLocations] = findpeaks(1-filt);

end

