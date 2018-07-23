function [ summary ] = summarizeMatrix( data, dimension, subsectionSize )
%SUMMARIZEMATRIX Summary of this function goes here
%   Detailed explanation goes here

if(~exist('subsectionSize', 'var'))
    subsectionSize = floor(size(data, dimension) / 4);
end

summary.meanValue = mean(data, dimension);
summary.stdValue = std(data, 0, dimension);
try
  summary.skewValue = skewness(data, 1, dimension);
catch ex
  summary.skewValue = ex;
end
try
  summary.kurtosisValue = kurtosis(data, 1, dimension);
catch ex
  summary.kurtosisValue = ex;
end
slice.startIndex = 1;
slice.endIndex = subsectionSize;
intervalCounter = 1;
while(slice.endIndex < size(data, dimension))
%   startIndex = min(find(cohInfo.timePoints >= slice.startTime));
%   endIndex = max(find(cohInfo.timePoints < slice.endTime));
  
  subsection = hslice(data, dimension, slice.startIndex:slice.endIndex);
  
  slice.meanValues = mean(subsection, dimension);
  slice.stdValues = std(subsection, 0, dimension);
  summary.slices(intervalCounter) = slice;
  intervalCounter = intervalCounter + 1;
  slice.startIndex = slice.startIndex + subsectionSize;
  slice.endIndex = slice.endIndex + subsectionSize;
end

end


