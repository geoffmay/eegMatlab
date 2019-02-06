function [ summary ] = asymptoteCoherenceReliability( icaCoh, frameCount )
%RESAMPLECOHERENCERELIABILITY Summary of this function goes here
%   Detailed explanation goes here

batchSize = 100;
queueLength = 10;
minDeltaFraction = 0.001;

totalFrames = size(icaCoh.matrix,1);


val1 = zeros(1, size(icaCoh.matrix, 2));
val2 = zeros(1, size(icaCoh.matrix, 2));
iterationCounter = 0;
delSum = zeros(1, size(icaCoh.matrix, 2));
delAvg = NaN(queueLength, size(icaCoh.matrix,2));
delCounter = 1;
finishedLoop = false;
deltaFractionBelowThreshold = false;
while(~deltaFractionBelowThreshold)
  %resample
  for repCounter = 1:batchSize;
    good = false;
    while(~good)
      start1 = ceil((totalFrames - frameCount) * rand);
      start2 = ceil((totalFrames - frameCount) * rand);
      %make sure samples don't overlap
      if(abs(start2 - start1) >= frameCount)
        try
          end1 = start1 + frameCount - 1;
          end2 = start2 + frameCount - 1;
          chunk1 = icaCoh.matrix(start1:end1, :);
          chunk2 = icaCoh.matrix(start2:end2, :);
          val1(1, :) = mean(chunk1, 1);
          val2(1, :) = mean(chunk2, 1);
          good = true;
        catch ex
        end
      end
    end
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
    fprintf('\niteration %d: avgDelFraction %f', iterationCounter, avgDelFraction);
    if(avgDelFraction < minDeltaFraction)
      deltaFractionBelowThreshold = true;
    end
  end
end
delStd = NaN(1,size(delArray,2));
checkDelAvg = NaN(1,size(delArray,2));
for i = 1:size(delArray,2)
  delStd(i) = std(delArray(:,i));
  checkDelAvg(i) = mean(delArray(:,i));
end

summary.frameCount = frameCount;
summary.iterations = iterationCounter;
summary.averageDifference = newAvg;
summary.averageDifferenceCheck = checkDelAvg;
summary.stddevDifference = delStd;
summary.finalDelDel = delDelFraction;
summary.queueLength  = queueLength;
summary.minDeltaFraction = minDeltaFraction;

end

