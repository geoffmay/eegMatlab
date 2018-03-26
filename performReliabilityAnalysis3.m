function [ summary ] = performReliabilityAnalysis3( eeg, outPath )
%PERFORMRELIABILITYANALYSIS Summary of this function goes here
%   Detailed explanation goes here



veryShort = 0; %flag for testing; loads a short file
makeDummy = 1; %creates a dummy file, which helps with running multiple instances
doMaster = 1; %combines all files
asymptote = 1; %flag to keep going until some type of convergence is reached.
injectNoise = 0;
doShuffle = 0;

%intermediateFile = sprintf('%s.intermediate.mat', outPath);

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

surfCoh = deriveRobiCoherenceMatrix(eeg);
if(injectNoise)
  eegNoise = eeg;
  %   means = nanmean(surfCoh.matrix, 1);
  scales = std(eeg.data, [], 2);
  injectionScale = 1;
  %   plot(means(1:end-1));
  for i = 1:length(scales)
    eegN = eegNoise.data(i, :);
    noise = wgn(size(eegN, 1), size(eegN,2), 0);
    eegNoise.data(i,:) = eegN + scales(i) .* noise;
  end
  noiseCoh = deriveRobiCoherenceMatrix(eegNoise);
end
if(doShuffle)
  permutation = randperm(size(surfCoh.matrix,1));
  permCoh.matrix = surfCoh.matrix(permutation, :);
  permCoh.labels = surfCoh.labels;
end
totalFrames = size(surfCoh.matrix,1);
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

%crunch the numbers (surface)
while(frameCount <= maxFrameCount)
  %   fprintf('\n(%s) resampling...',char(datetime));
  %
  %
  %   fprintf('\n%s resampling segments of length %d', char(datetime), frameCount);
  if(asymptote)
    %     summary.icaResample(lengthCounter) = asymptoteCoherenceReliability2(icaCoh, frameCount);
    summary.surfResample(lengthCounter) = asymptoteCoherenceReliability2(surfCoh, frameCount);
    if(injectNoise)
      summary.surfNoiseResample(lengthCounter) = asymptoteCoherenceReliability2(noiseCoh, frameCount);
    end
    if(doShuffle)
      summary.surfPermResample(lengthCounter) = asymptoteCoherenceReliability2(permCoh, frameCount);
    end
    %     save(intermediateFile, '-v7.3');
  end
  
  frameCount = frameCount + frameSkip;
  lengthCounter = lengthCounter + 1;
  save(intermediateFile);
  fprintf('...done');
end
% % delete(intermediateFile);
% summary.icaLabels = icaCoh.labels;
summary.surfLabels = surfCoh.labels;
% summary.icaInfo = icaCoh.icaInfo;

clear surfCoh;
clear noiseCoh;
clear permCoh;

icaCoh = deriveIcaCoherenceMatrix(eeg, 128);
% totalFrames = size(icaCoh.matrix,1);

if(injectNoise)
%   eegNoise = eeg;
%   %   means = nanmean(surfCoh.matrix, 1);
%   scales = std(eeg.data, [], 2);
%   injectionScale = 1;
%   %   plot(means(1:end-1));
%   for i = 1:length(scales)
%     eegN = eegNoise.data(i, :);
%     noise = wgn(size(eegN, 1), size(eegN,2), 0);
%     eegNoise.data(i,:) = eegN + scales(i) .* noise;
%   end
  noiseCoh = deriveIcaCoherenceMatrix(eegNoise);
end
if(doShuffle)
%   permutation = randperm(size(surfCoh.matrix,1));
  permCoh.matrix = icaCoh.matrix(permutation, :);
  permCoh.labels = icaCoh.labels;
  
end

lengthCounter = 1;
frameCount = minFrameCount;

%crunch the numbers (ica)
while(frameCount <= maxFrameCount)
  fprintf('\n(%s) resampling...',char(datetime));
  
  
  fprintf('\n%s resampling segments of length %d', char(datetime), frameCount);
  if(asymptote)
    summary.icaResample(lengthCounter) = asymptoteCoherenceReliability2(icaCoh, frameCount);
    if(injectNoise)
      summary.icaNoiseResample(lengthCounter) = asymptoteCoherenceReliability2(noiseCoh, frameCount);
    end
    if(doShuffle)
      summary.icaPermResample(lengthCounter) = asymptoteCoherenceReliability2(permCoh, frameCount);
    end
    %     summary.surfResample(lengthCounter) = asymptoteCoherenceReliability2(surfCoh, frameCount);
    %     save(intermediateFile, '-v7.3');
  end
  
  frameCount = frameCount + frameSkip;
  lengthCounter = lengthCounter + 1;
  save(intermediateFile);
  fprintf('...done');
end
% delete(intermediateFile);
summary.icaLabels = icaCoh.labels;
summary.icaInfo = icaCoh.icaInfo;

clear icaCoh;

%save('/home/data/EEG/processed/Oregon/workspace.mat');
