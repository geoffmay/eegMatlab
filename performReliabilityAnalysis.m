function [ summary ] = performReliabilityAnalysis( eeg, outPath )
%PERFORMRELIABILITYANALYSIS Summary of this function goes here
%   Detailed explanation goes here



veryShort = 0; %flag for testing; loads a short file
makeDummy = 1; %creates a dummy file, which helps with running multiple instances
doMaster = 1; %combines all files
if(~veryShort)
  repetitions = 1000;
  minFrameCount = 128 * 15;
  frameSkip = 128 * 15;
  maxFrameCount = 128 * 60 * 25;
else
  repetitions = 10;
  minFrameCount = 10;
  frameSkip = 10;
  maxFrameCount = 100;
  
end

icaCoh = deriveIcaCoherenceMatrix(eeg, 128);
totalFrames = size(icaCoh.matrix,1);
if(totalFrames < 128 * 120) %less than two minutes
  minFrameCount = ceil(totalFrames / 100);
  frameSkip = minFrameCount;
end
maxFrameCount = min(maxFrameCount, floor(totalFrames * .4));
%iterate through frame counts to give greater and greater lengths
frameCount = minFrameCount;
icaRhos = NaN(floor(maxFrameCount / frameSkip), size(icaCoh.matrix, 2));
lengthCounter = 1;
while(frameCount <= maxFrameCount)
  summary.sampleFrameDuration(lengthCounter) = frameCount;
  fprintf('\n%s resampling segments of length %d', char(datetime), frameCount);
  val1 = NaN(repetitions, size(icaCoh.matrix, 2));
  val2 = NaN(repetitions, size(icaCoh.matrix, 2));
  %resample
  for repCounter = 1:repetitions
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
          val1(repCounter, :) = mean(chunk1, 1);
          val2(repCounter, :) = mean(chunk2, 1);
          good = true;
        catch ex
        end
      end
    end
  end
  %correlate
  fprintf('...%s correlating', char(datetime));
  for varCounter = 1:size(val1,2)
    x = val1(:, varCounter);
    y = val2(:, varCounter);
    icaRhos(lengthCounter, varCounter) = corr(x, y);
  end
  frameCount = frameCount + frameSkip;
  lengthCounter = lengthCounter + 1;
  fprintf('...done');
end

summary.icaRhos = icaRhos;
summary.icaLabels = icaCoh.labels;
summary.icaWeightInfo = icaCoh.icaInfo;


surfInfo = deriveRobiCoherenceMatrix(eeg);

%again for surface coherecne, iterate through frame counts
frameCount = minFrameCount;
surfRhos = NaN(floor(maxFrameCount / frameSkip), size(icaCoh.matrix, 2));
% surfRhos = NaN(maxFrameCount / frameSkip, size(icaCoh.matrix, 2));
lengthCounter = 1;
while(frameCount <= maxFrameCount)
  fprintf('\n%s resampling segments of length %d', char(datetime), frameCount);
  val1 = NaN(repetitions, size(surfInfo.matrix, 2));
  val2 = NaN(repetitions, size(surfInfo.matrix, 2));
  %resample
  for repCounter = 1:repetitions
    good = false;
    while(~good)
      start1 = ceil((totalFrames - frameCount) * rand);
      start2 = ceil((totalFrames - frameCount) * rand);
      %make sure samples don't overlap
      if(abs(start2 - start1) >= frameCount)
        try
          end1 = start1 + frameCount - 1;
          end2 = start2 + frameCount - 1;
          chunk1 = surfInfo.matrix(start1:end1, :);
          chunk2 = surfInfo.matrix(start2:end2, :);
          val1(repCounter, :) = mean(chunk1, 1);
          val2(repCounter, :) = mean(chunk2, 1);
          good = true;
        catch ex
        end
      end
    end
  end
  %correlate
  for varCounter = 1:size(val1,2)
    x = val1(:, varCounter);
    y = val2(:, varCounter);
    surfRhos(lengthCounter, varCounter) = corr(x, y);
  end
  frameCount = frameCount + frameSkip;
  lengthCounter = lengthCounter + 1;
end

summary.surfRhos = surfRhos;
summary.surfLabels = surfInfo.labels;
summary.meanSurfRho = mean(surfRhos, 2);
summary.meanIcaRho = mean(icaRhos, 2);
if(false)
  figure;
  plot([summary.meanSurfRho, summary.meanIcaRho]);
  legend({'surface', 'ica'});
  
  scatter(x,y);
end

end

