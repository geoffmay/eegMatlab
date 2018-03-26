 
function [ summary ] = asymptoteCoherenceReliability2( surfCoh, frameCount )
%RESAMPLECOHERENCERELIABILITY Summary of this function goes here
%   Detailed explanation goes here

batchSize = 100;
queueLength = 10;
minDeltaFraction = 0.001;

%(31-Jan-2018 06:10:36) iteration 1000: avgDelFraction 0.015277
%(31-Jan-2018 06:12:25) iteration 7600: avgDelFraction 0.000992
%(31-Jan-2018 07:26:49) iteration 56000: avgDelFraction 0.000100

totalFrames = size(surfCoh.matrix,1);

iterationCounter = 0;
delSum = zeros(1, size(surfCoh.matrix, 2));
delAvg = NaN(queueLength, size(surfCoh.matrix,2));
delCounter = 1;
finishedLoop = false;
deltaFractionBelowThreshold = false;

fprintf('\n(%s): computing subsample averages length %d', char(datetime), frameCount);
%(new & fast): compute all averages beforehand
starts = 1:totalFrames-frameCount + 1;
ends = frameCount:totalFrames;
avgSurfMat = NaN(length(starts), size(surfCoh.matrix,2));
val1 = zeros(1, size(avgSurfMat, 2));
val2 = zeros(1, size(avgSurfMat, 2));

total = zeros(1, size(surfCoh.matrix,2));
recip = 1/frameCount;
for i = 1:totalFrames
  i2 = i - frameCount+1;
  total = total + surfCoh.matrix(i, :);
  if(i2 > 0)
    total = total - surfCoh.matrix(i2, :);
    avgSurfMat(i2,:) = total .* recip;
  end
end

%fprintf('\n(%s): randomly resampling subsaples', char(datetime));
%      %to choose subsampel location intelligently, make sure that a second
      %contiguous subsample exists on one side or the other.
maxTroubleFreeSubsampleLength = floor((totalFrames+2)/3);
start1Mask = logical(ones(1, length(starts)));
if(frameCount > maxTroubleFreeSubsampleLength)   
  hasLeft = starts > frameCount;
  hasRight = starts < (length(starts)-frameCount);
  start1Mask = start1Mask & (hasLeft | hasRight);
end
start1Ind = find(start1Mask);

while(~deltaFractionBelowThreshold)
  %resample
  for repCounter = 1:batchSize;
    %     good = false;
    %     while(~good)
    start1a = ceil(length(start1Ind) * rand);
    start1 = start1Ind(start1a);
    start2Mask = start1Mask;
    invalidStart2a = max(1, start1 - frameCount + 1);
    invalidStart2b = min(length(starts), start1 + frameCount - 1);
    start2Mask(invalidStart2a:invalidStart2b) = 0;
    start2Ind = find(start2Mask);
    start2a = ceil(length(start2Ind) * rand);
    start2 = start2Ind(start2a);
    
    start2 = ceil((totalFrames - frameCount) * rand);
    %       %make sure samples don't overlap
    %       if(abs(start2 - start1) >= frameCount)
    %         try
    %           end1 = start1 + frameCount - 1;
    %           end2 = start2 + frameCount - 1;
    %           chunk1 = surfCoh.matrix(start1:end1, :);
    %           chunk2 = surfCoh.matrix(start2:end2, :);
    %       val1(1, :) = mean(chunk1, 1);
    %       val2(1, :) = mean(chunk2, 1);
    %           good = true;
    %         catch ex
    %         end
    %     end
    val1 = avgSurfMat(start1,:);
    val2 = avgSurfMat(start2,:);
    
    %     end
    iterationCounter = iterationCounter + 1;
    del = abs(val1 - val2);
    delArray(iterationCounter, :) = del;
    delSum = delSum + del;
  end
  newAvg = delSum ./ iterationCounter;
  delAvg(delCounter, :) = newAvg;
  delCounter = delCounter + 1;
  if(delCounter > queueLength)
    delCounter = 1;
    finishedLoop = true;
  end
  if(finishedLoop)
    delDel = abs(diff(delAvg, 1, 1));
    meanDelDel = mean(delDel, 1);
    meanDel = mean(delAvg);
    delDelFraction = meanDelDel ./ meanDel;
    avgDelFraction = mean(delDelFraction, 2);
%     fprintf('\n(%s) iteration %d: avgDelFraction %f', char(datetime), iterationCounter, avgDelFraction);
    if(avgDelFraction < minDeltaFraction)
      deltaFractionBelowThreshold = true;
    end
  end
end
delStd = NaN(1,size(delArray,2));
checkDelAvg = NaN(1,size(delArray,2));
prc.percentileKeys = [1, sqrt(.05)*100, 5, 95, sqrt(.95)*100, 99];
prc.percentileValues = NaN(length(prc.percentileKeys), size(delArray,2));
for i = 1:size(delArray,2)
  delStd(i) = std(delArray(:,i));
  checkDelAvg(i) = mean(delArray(:,i));
  prc.percentileValues(:,i) = prctile(delArray(:,i), prc.percentileKeys);
end

summary.frameCount = frameCount;
summary.iterations = iterationCounter;
summary.averageDifference = newAvg;
summary.averageDifferenceCheck = checkDelAvg;
summary.stddevDifference = delStd;
summary.finalDelDel = delDelFraction;
summary.queueLength  = queueLength;
summary.minDeltaFraction = minDeltaFraction;
summary.percentiles = prc;

end

