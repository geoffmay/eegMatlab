function [ summary ] = performReliabilityAnalysis1( eeg, outPath )
%PERFORMRELIABILITYANALYSIS Summary of this function goes here
%   Detailed explanation goes here



veryShort = 0; %flag for testing; loads a short file
makeDummy = 1; %creates a dummy file, which helps with running multiple instances
doMaster = 1; %combines all files
asymptote = 1; %flag to keep going until some type of convergence is reached.

intermediateFile = sprintf('%s.intermediate.mat', outPath);

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
  eeg.data(:, 10000:end) = [];
  
end

icaCoh = deriveIcaCoherenceMatrix(eeg, 128);
surfCoh = deriveRobiCoherenceMatrix(eeg);

totalFrames = size(icaCoh.matrix,1);
if(totalFrames < 128 * 120) %less than two minutes
  minFrameCount = ceil(totalFrames / 100);
  frameSkip = minFrameCount;
end
maxFrameCount = min(maxFrameCount, floor(totalFrames * .4));
%iterate through frame counts to give greater and greater lengths
% icaRhos = NaN(floor(maxFrameCount / frameSkip), size(icaCoh.matrix, 2));
frameCount = minFrameCount;
lengthCounter = 0;
while(frameCount <= maxFrameCount)
  frameCount = frameCount + frameSkip;
  lengthCounter = lengthCounter + 1;
end
lengthCounter = 1;
frameCount = minFrameCount;

%crunch the numbers
while(frameCount <= maxFrameCount)
  fprintf('\n%s resampling segments of length %d', char(datetime), frameCount);
  if(asymptote)
    summary.icaResample(lengthCounter) = asymptoteCoherenceReliability(icaCoh, frameCount);
    summary.surfResample(lengthCounter) = asymptoteCoherenceReliability(surfCoh, frameCount);
    save(intermediateFile, '-v7.3');
  end
  
  frameCount = frameCount + frameSkip;
  lengthCounter = lengthCounter + 1;
  fprintf('...done');
end
delete(intermediateFile);
summary.icaLabels = icaCoh.labels;
summary.surfLabels = surfCoh.labels;

%salvage
if(false)
  intermediateFile = sprintf('%s.intermediate.mat', outPath);
  intermediate = load(intermediateFile)
  intermediate.summary
  intermediate.summary.icaLabels = intermediate.icaCoh.labels;
  intermediate.summary.surfLabels = intermediate.surfCoh.labels;
  summary = intermediate.summary;
  save(outPath, 'summary');
end


%analysis section
if(false)
  for i = 1:length(summary.icaResample)
    icaDiffs(i,:) = summary.icaResample(i).averageDifference;
  end
  for i = 1:length(summary.surfResample)
    surfDiffs(i,:) = summary.surfResample(i).averageDifference;
  end
  meanSurfDiff = mean(surfDiffs, 2);
  meanIcaDiff = mean(icaDiffs, 2);
  figure;
  hold on;
  ind = 1:10;
  x = summary.sampleFrameDuration ./ (128 * 60); x = x';
  plot(x, surfDiffs(:,ind));
  legend(summary.surfLabels(ind));
end

end

