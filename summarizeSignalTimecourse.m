function [summary ] = summarizeSignalTimecourse(cohInfo, intervalSize)
%SUMMARIZESIGNALTIMECOURSE Summary of this function goes here
%   Detailed explanation goes here

if(~exist('intervalSize', 'var'))
    intervalSize = 128 * 60 * 0.25;
end

summary.meanValue = mean(cohInfo);
summary.stdValue = std(cohInfo);
summary.skewValue = skewness(cohInfo);
summary.kurtosisValue = kurtosis(cohInfo);
slice.startFrame = 1;
slice.endFrame = intervalSize;
intervalCounter = 1;
while(slice.endFrame < length(cohInfo))
%   startIndex = min(find(cohInfo.timePoints >= slice.startFrame));
%   endIndex = max(find(cohInfo.timePoints < slice.endFrame));
  data = cohInfo(slice.startFrame:slice.endFrame);
  slice.meanValues = mean(data);
  slice.stdValues = std(data);
  summary.slices(intervalCounter) = slice;
  intervalCounter = intervalCounter + 1;
  slice.startFrame = slice.startFrame + intervalSize;
  slice.endFrame = slice.endFrame + intervalSize;
end

end


