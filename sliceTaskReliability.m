% inputFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\Oregon\artifactRemovedAnalytics';
% oldFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\Oregon\artifactRemovedReliability';

fewMeasures = false;
includeOldData = false;
doPlot = false;
inputFolder = '/home/data/EEG/processed/Oregon/artifactRemoved';
outputFolder = '/home/data/EEG/processed/Oregon/separatedReliability';



inputFiles = dir(inputFolder);
inputFiles([inputFiles.isdir]) = [];

% inputFile = 'PM101Surface.mat';


for fileCounter = 1:length(inputFiles)
  inputFile = inputFiles(fileCounter).name;
  outputPath = fullfile(outputFolder, inputFile);
  
  if(~exist(outputPath, 'var'))
    load(fullfile(inputFolder, inputFile));
    
    task = guessWahbehTaskFromAlpha(surfCoh);
    
    eyesClosed = surfCoh.matrix(task == 1, :);
    eyesOpen = surfCoh.matrix(task == -1, :);
    
    %compute reliability for eyes open and eyes closed states
    minDeltaFraction = 1e-3;
    i = 1;
    frameStep = 5;
    frameCount = 128*frameStep*i;
    while(frameCount < (length(eyesOpen)*.4) || frameCount < (length(eyesClosed)*.4))
      reliabilities.frameCounts(i) = frameCount;
      if( frameCount < (length(eyesClosed)*.4))
        reliabilities.eyesClosed(i) = asymptoteSignalReliability(eyesClosed, frameCount);
        fprintf(' (eyes closed)')
      end
      if( frameCount < (length(eyesOpen)*.4))
        reliabilities.eyesOpen(i) = asymptoteSignalReliability(eyesOpen, frameCount);
        fprintf(' (eyes open)')
      end
      i = i + 1;
      frameCount = 128*5*i;
    end
    
    if(includeOldData)
      %load old combined reliability data
      oldData = load(fullfile(oldFolder, inputFile));
      oldInd = find(strcmp(oldData.summary.surfLabels, targetMeasure));
      stdCoeff = norminv(.95);
      for i = 1:length(oldData.summary.surfResample)
        a = oldData.summary.surfResample(i);
        oldPlot(i) = a.averageDifference(oldInd) + stdCoeff * a.stddevDifference(oldInd);
        xOld(i) = a.frameCount / 128;
      end
      
      a = [reliabilities.eyesOpen];
      a = [a.percentiles];
      for i = 1:length(a)
        eyesOpen95(i) = a(i).percentileValues(3);
      end
      xOpen = (1:length(eyesOpen95)) .* frameStep;
      
      a = [reliabilities.eyesClosed];
      a = [a.percentiles];
      for i = 1:length(a)
        eyesClosed95(i) = a(i).percentileValues(3);
      end
      xClosed = (1:length(eyesClosed95)) .* frameStep;
    end
    
    
    %plot reliabilities
    if(doPlot)
      figure;
      hold on;
      for i = 1:length(plotX)
        plot(plotX{i}, plotY{i});
      end
      legend(plotLabel);
      xlabel('subsample duration (seconds)');
      ylabel('power difference (log(uv^2))');
      %     plot(xClosed, eyesClosed95);
      %     plot(xOpen, eyesOpen95);
      %     plot(xOld, oldPlot);
      %     legend({'alpha closed','alpha open','alpha combined'});
    end
    
    
    % save('C:\Users\Neuro\Documents\MATLAB\processed\Oregon\ptsdTaskSliceDemo.mat');
    
    save(outputPath, 'reliabilities');
    clear reliabilities
  end
end
