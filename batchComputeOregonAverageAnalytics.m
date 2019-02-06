inputFolder = '/home/data/EEG/data/Oregon/PTSD';
outputFolder = '/home/data/EEG/processed/Oregon/avgValues';


veryShort = 0; %flag for testing; loads a short file
makeDummy = 1; %creates a dummy file, which helps with running multiple instances
doMaster = 1; %combines all files
repetitions = 1000;
minFrameCount = 128 * 15;
frameSkip = 128 * 15;
maxFrameCount = 128 * 60 * 25;

fileNames = dir(inputFolder);
fileNames([fileNames.isdir]) = [];

for fileCounter = 1:length(fileNames);
  fileName = fullfile(inputFolder, fileNames(fileCounter).name);
  [folder, file, ext] = fileparts(fileName);
  outFile = [file '.mat'];
  outPath = fullfile(outputFolder, outFile);
  
  eeg = loadBdf(fileName);
  processThisFile = true;
  if(makeDummy)
    if(~exist(outPath, 'file'))
      placeHolder = sprintf('analysis started on %s', char(datetime));
      save(outPath, 'placeHolder');
      processThisFile = true;
    else
      processThisFile = false;
    end
  end
  if(processThisFile)
    
    icaCoh = deriveIcaCoherenceMatrix(eeg, 128);
    surfCoh = deriveRobiCoherenceMatrix(eeg);
    
    surfSummary.meanValue = mean(surfCoh.matrix,1);
    surfSummary.stdValue = std(surfCoh.matrix,1);
    surfSummary.skewValue = skewness(surfCoh.matrix,1);
    surfSummary.kurtosisValue = kurtosis(surfCoh.matrix,1);
    surfSummary.labels = surfCoh.labels;    
    icaSummary.meanValue = mean(icaCoh.matrix,1);
    icaSummary.stdValue = std(icaCoh.matrix,1);
    icaSummary.skewValue = skewness(icaCoh.matrix,1);
    icaSummary.kurtosisValue = kurtosis(icaCoh.matrix,1);
    icaSummary.labels = icaCoh.labels;
    icaSummary.icaInfo = icaCoh.icaInfo;    

    icaSummary.note = 'time is denoted in minutes';
    surfSummary.note = 'time is denoted in minutes';
    intervalSize = 0.25;

    icaSlice.startTime = 0;
    icaSlice.endTime = intervalSize;
    surfSlice.startTime = 0;
    surfSlice.endTime = intervalSize;
    intervalCounter = 1;
    while(icaSlice.endTime < max(icaCoh.timePoints))
      startIndex = min(find(icaCoh.timePoints >= icaSlice.startTime));
      endIndex = max(find(icaCoh.timePoints < icaSlice.endTime));
      icaSlice.meanValues = mean(icaCoh.matrix(startIndex:endIndex,:, 1));
      icaSlice.stdValues = std(icaCoh.matrix(startIndex:endIndex,:, 1));
      surfSlice.meanValues = mean(surfCoh.matrix(startIndex:endIndex,:, 1));
      surfSlice.stdValues = std(surfCoh.matrix(startIndex:endIndex,:, 1));
      icaSummary.slices(intervalCounter) = icaSlice;
      surfSummary.slices(intervalCounter) = surfSlice;
      intervalCounter = intervalCounter + 1;
      icaSlice.startTime = icaSlice.startTime + intervalSize;
      icaSlice.endTime = icaSlice.endTime + intervalSize;      
      surfSlice.startTime = surfSlice.startTime + intervalSize;
      surfSlice.endTime = surfSlice.endTime + intervalSize;      
    end
    
    summary.icaSummary = icaSummary;
    summary.surfSummary = surfSummary;
    save(outPath, 'summary');
    
    if(false)
      slices = summary.icaSummary.slices;
      means = NaN(length(slices), length(slices(1).meanValues));
      stds = NaN(length(slices), length(slices(1).meanValues));
      for i = 1:length(slices)
        means(i,:) = slices(i).meanValues;
        stds(i,:) = slices(i).stdValues;
      end
    end
  end
end