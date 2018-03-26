function [ badReferenceCutoff ] = checkForBadReference( coherenceTimeCourse )
%CHECKFORBADREFERENCE tells the sample index at which coherence plateaus
%around 1, which indicates that a reference (such as mastoid) has fallen
%off.
%   This assumes a sample rate of 128.

coh = coherenceTimeCourse;
windowSize = 1024;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
% try
coh1 = filtfilt(b, a, coh);
threshold1 = 0.99;
supraThreshold = coh1 > threshold1;
supraThresholdIndexes = find(supraThreshold);
if(length(supraThresholdIndexes) == 0)
    badReferenceCutoff = length(coh1);
else
    threshold2 = 0.9;
    minPercentBelowThreshold2 = .1;
    
    thresholdMetaIndex = 0;
    while(percentBelowThreshold < minPercentBelowThreshold & thresholdMetaIndex <= length(supraThresholdIndexes))
        thresholdMetaIndex = thresholdMetaIndex + 1;
        cutoffPoint = supraThresholdIndexes(thresholdMetaIndex);
        subThreshold = coh(cutoffPoint:end) < threshold2;
        percentBelowThreshold = sum(subThreshold) / (length(coh) - cutoffPoint + 1);
    end
    subThreshold = coh < threshold2;
    subThreshold(cutoffPoint:end) = [];
    lastUnfilteredSubThreshold = find(subThreshold);
    if(length(lastUnfilteredSubThreshold) > 0)
        lastUnfilteredSubThreshold = lastUnfilteredSubThreshold(end);
        badReferenceCutoff = lastUnfilteredSubThreshold;
    else
        badReferenceCutoff = 1;
    end
end

% catch err
%     badReferenceCutoff = 1;
% end

end

