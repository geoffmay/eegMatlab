coherenceFilenames = [...
  {'/home/data/EEG/processed/Robi/coherence/ROBI_003_baseline eyes open_630158995692243270PlusPhase.mat'},...
  {'/home/data/EEG/processed/Robi/coherence/ROBI_003_outcome eyes open_630230539149591228PlusPhase.mat'}...
  ];

outputFilenames = [...
  {'/home/data/EEG/processed/Robi/coherence/ROBI_003_baseline eyes open_highCoherencePhaseSlope.mat'},...
  {'/home/data/EEG/processed/Robi/coherence/ROBI_003_outcome eyes open_highCoherencePhaseSlope.mat'}...
  ];

outputFolder = '/home/data/EEG/processed/Robi/coherence/highCoherencePhaseLag/';


for fileCounter = 1:length(coherenceFilenames);
  coherenceFilename = coherenceFilenames{fileCounter};
  outputPath = outputFilenames{fileCounter};
  if(~exist(outputPath))
    load(coherenceFilename);
    outputPath = getRobiOutputFilenames({filename}, outputFolder);
    outputPath = outputPath{1};
    
    fileInfo = dir(filename);
    fileId = fopen(filename);
    contents = fread(fileId, fileInfo.bytes/8, 'double');
    fclose(fileId);
    channelCount = 34;
    fileLength = fileInfo.bytes/8/channelCount;
    data = reshape(contents, 34, fileLength)';
    clear contents;
    
    pair = channelPairs(1);
    deltaThreshold = 0.9;
    minEpochLength = 129;
    
    frequencySelector = [1 4];
    normalizedPhaseAngle = pair.phaseAngle(:,frequencySelector) .* (1/(2*pi)) + 0.5;
    plotMatrix = [pair.coherence(:,frequencySelector), normalizedPhaseAngle];
    x = (1:size(plotMatrix,1)) ./ 2048;
    
    deltaCoh = NaN(size(channelPairs(1).coherence, 1), length(channelPairs));
    
    for fileCounter = 1:size(deltaCoh,2)
      a = [channelPairs(fileCounter).coherence];
      deltaCoh(:,fileCounter) = a(:,1);
    end
    supraThreshold = deltaCoh > deltaThreshold;
    for fileCounter = 1:size(deltaCoh,2)
      fprintf('.');
      if(mod(fileCounter, 100)==0)
        fprintf('\n%d', fileCounter);
      end
      a = supraThreshold(:,fileCounter);
      b = diff(a);
      rising = find(b == 1);
      falling = find(b == -1);
      while(falling(1) < rising(1))
        falling(1) = [];
      end
      while(length(rising) > length(falling))
        rising(end) = [];
      end
      while(length(rising) < length(falling))
        falling(end) = [];
      end
      %chop the end off because coherence has some lag in finding quality
      %data.
      falling1 = falling - 64;
      falling = falling1;
      epochLengths = falling - rising;
      goodEpoch = find(epochLengths > minEpochLength);
      clear slopes phasePlots polyfits;
      [chan1, chan2] = channelsFromPairIndex(fileCounter, size(deltaCoh,2));
      for j = 1:length(goodEpoch)
        startIndex = rising(goodEpoch(j));
        endIndex = falling(goodEpoch(j));
        dataStartIndex = startIndex * 16 + 1024;
        dataEndIndex = endIndex * 16 + 1024;
        [chanPair.slopes(j), chanPair.phasePlots{j}, chanPair.polyfits{j}] = phaseSlopeIndex(data(dataStartIndex:dataEndIndex,:), chan1, chan2);
      end
      if(false)
        close all;
        hist(slopes);
        close all;
        figure;
        hold on;
        for j = 1:length(phasePlots)
          plot(phasePlots{j});
        end
      end
      chanPairs(fileCounter) = chanPair;
    end
    save(outputPath, 'filename', 'coherenceFilename', 'chanPairs', '-v7.3');
  end
end

%plot coherence;
if(false)
  figure;
  plot(x, plotMatrix);
  if(length(frequencySelector) <= 1)
    title(sprintf('%s %s', pair.label, freqLabels{frequencySelector}));
    legend('coherence', 'phaseAngle');
  else
    title(sprintf('%s', pair.label, freqLabels{frequencySelector}));
    legendCounter = 1;
    legendText = cell(0);
    for i = 1:length(frequencySelector)
      freqText = freqLabels{frequencySelector(i)};
      legendText{legendCounter} = sprintf('%s coherence', freqText);
      legendCounter = legendCounter + 1;
    end
    for i = 1:length(frequencySelector)
      freqText = freqLabels{frequencySelector(i)};
      legendText{legendCounter} = sprintf('%s phase angle', freqText);
      legendCounter = legendCounter + 1;
    end
    legend(legendText);
    
  end
  zoom xon;
  pan xon;
  
end
% else
%   if(~exist('baseline', 'var'))
%     baseline = load(outputFilenames{1});
%     outcome = load(outputFilenames{2});
%   end
%     labels = antChannelLocs;
%   for pairCount = 1:496
%     [i, j] = channelsFromPairIndex(pairCount, 496);
%     label = sprintf('%s-%s', labels{i}, labels{j});
%     close all;
%     figure;
%     hist([baseline.chanPairs.slopes]);
%     title(sprintf('%s baseline', label));
%     figure;
%     hist([outcome.chanPairs.slopes]);
%     title(sprintf('%s outcome', label));
%     tilefigs;
%   end
% end

for fileCounter = 1:length(outputFilenames)
  data = load(outputFilenames{fileCounter});
  for i = 1:length(data.chanPairs)
    for j = 1:length(data.chanPairs(i).phasePlots)
      phasePlot = data.chanPairs(i).phasePlots{j};
      if(j==1)
        phasePlotSum = zeros(size(phasePlot));
      end
      if(all(~isnan(phasePlot)))
        phasePlotSum = phasePlotSum + phasePlot;        
      end
    end
    close all;
    plot(phasePlotSum);
  end
end
