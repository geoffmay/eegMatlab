function [ summary ] = summarizeCoherenceMatrix( cohInfo, intervalSize )
%SUMMARIZECOHERENCE Summary of this function goes here
%   Detailed explanation goes here

if(~exist('intervalSize', 'var'))
    intervalSize = 0.25;
end

summary.meanValue = mean(cohInfo.matrix,1);
summary.stdValue = std(cohInfo.matrix,1);
summary.skewValue = skewness(cohInfo.matrix,1);
summary.kurtosisValue = kurtosis(cohInfo.matrix,1);
summary.labels = cohInfo.labels;
slice.startTime = 0;
slice.endTime = intervalSize;
intervalCounter = 1;
while(slice.endTime < max(cohInfo.timePoints))
  startIndex = min(find(cohInfo.timePoints >= slice.startTime));
  endIndex = max(find(cohInfo.timePoints < slice.endTime));
  slice.meanValues = mean(cohInfo.matrix(startIndex:endIndex,:, 1));
  slice.stdValues = std(cohInfo.matrix(startIndex:endIndex,:, 1));
  summary.slices(intervalCounter) = slice;
  intervalCounter = intervalCounter + 1;
  slice.startTime = slice.startTime + intervalSize;
  slice.endTime = slice.endTime + intervalSize;
end

end

