function [ badReferenceCutoff ] = checkForBadReference( coherenceTimeCourse )
%CHECKFORBADREFERENCE tells the sample index at which coherence plateaus
%around 1, which indicates that a reference (such as mastoid) has fallen
%off.
%   This assumes a sample rate of 128.

coh = coherenceTimeCourse;
windowSize = 1024;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
try
    coh1 = filtfilt(b, a, coh);
catch err
    badReferenceCutoff = 1;
end
threshold1 = 0.99;
supraThreshold = coh1 > threshold1;
firstFilteredSupraThreshold = find(supraThreshold);
if(length(firstFilteredSupraThreshold) == 0)
    badReferenceCutoff = length(coh1);
else
    firstFilteredSupraThreshold = firstFilteredSupraThreshold(1);
    threshold2 = 0.9;
    subThreshold = coh < threshold2;
    subThreshold(firstFilteredSupraThreshold:end) = [];
    lastUnfilteredSubThreshold = find(subThreshold);
    lastUnfilteredSubThreshold = lastUnfilteredSubThreshold(end);
    badReferenceCutoff = lastUnfilteredSubThreshold;
end


end

