clear;
verbose = false;
unsortedFiles = getRobiDataFiles;
eyesOpen = cellfun(@length, strfind(unsortedFiles, 'eyes open')) > 0;
files = [unsortedFiles(eyesOpen), unsortedFiles(~eyesOpen)];
windowSize = 10000;

for fileCounter = 1:length(files)
  %fileCounter = 1
  filename = files{fileCounter};
  
  outputFolder = '/home/data/EEG/processed/Robi/summaries';
  outputFilename = strrep(filename, '/home/data/EEG/data/ROBI/', '');
  outputFilename = strrep(outputFilename, '/', '.');
  outputFilename = sprintf('%s-%d.mat', outputFilename, windowSize);
  outputPath = fullfile(outputFolder, outputFilename);
  
  if(~exist(outputPath, 'file'))
    placeholder = sprintf('started on %s', char(datetime));
    save(outputPath, 'placeholder');
    clear('phaseSlopes');
    
    %phase slope topography
    %look at different slices
    eeg = loadRobiEeg(filename);
    oldEeg = eeg;
    counter = 1;
    counterMax = floor(size(oldEeg.data,2) / windowSize)-1;
    dataStart = windowSize*counter;
    dataEnd = windowSize*(counter+1);
    while(dataEnd < size(oldEeg.data,2))
      fprintf('\n%s: %f', char(datetime), dataEnd / size(oldEeg.data,2) * 100);
      eeg.data = oldEeg.data(:, dataStart:dataEnd);
      phaseSlope = phaseSlopeTopography(eeg);
      slice.lags = phaseSlope.estimatedTimeLag;
      slice.snr = phaseSlope.signalToNoiseRatio;
      slice.startTime = eeg.times(dataStart);
      slice.endTimes = eeg.times(dataEnd);
      phaseSlopes(counter) = slice;
      counter = counter + 1;
      dataStart = windowSize*counter;
      dataEnd = windowSize*(counter+1);
    end
    
    verbose = false;
    if(verbose)
      
      topoplot(lags(1,:), eeg.chanlocs);
      
      figure;
      imagesc(lags')
      colorbar
      set(gca, 'ytick', 1:eeg.nbchan);
      ytickLabels = {eeg.chanlocs.labels};
      set(gca, 'ytickLabel', ytickLabels);
      
      figure;
      plot(lags);
      
    end
    
    [coherenceData, measureSummary] = deriveRobiCoherenceMatrix(filename);
    try
      save(outputPath, 'phaseSlopes', 'measureSummary', '-v7.3');
    catch ex
      %do nothing
    end
    
  end
end
