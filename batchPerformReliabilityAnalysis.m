
finishIncompleteFiles = 0;
useNewAsymptotes = 1;

veryShort = 0; %flag for testing; loads a short file
makeDummy = 1; %creates a dummy file, which helps with running multiple instances
doMaster = 1; %combines all files
repetitions = 1000;
minFrameCount = 128 * 15;
frameSkip = 128 * 15;
maxFrameCount = 128 * 60 * 25;



inputFolder = '/home/data/EEG/data/Oregon/PTSD';
inputFolder = '/home/data/EEG/processed/Oregon/coherence';
outputFolder = '/home/data/EEG/processed/Oregon/reliability8';
outputFolder = '/home/data/EEG/processed/Oregon/reliability5';
outputFolder = '/home/data/EEG/processed/Oregon/reliability9';
if(~exist(outputFolder))
  mkdir(outputFolder);
end
fileNames = dir(inputFolder);
fileNames([fileNames.isdir]) = [];


%loop through all files
for fileCounter = 1:length(fileNames)
  fileName = fullfile(inputFolder, fileNames(fileCounter).name);
  [folder, file, ext] = fileparts(fileName);
  outFile = [file '.mat'];
  outPath = fullfile(outputFolder, outFile);
  
  processThisFile = true;
  if(makeDummy)
    if(~exist(outPath, 'file'))
      placeHolder = sprintf('analysis started on %s', char(datetime));
      save(outPath, 'placeHolder');
      processThisFile = true;
    else
      if(~finishIncompleteFiles)
        processThisFile = false;
      else
        [tempFolder, tempFile, tempExt] = fileparts(fileName);
        intermediateFilename = fullfile(outputFolder, [tempFile '.mat.intermediate.mat']);
        if(exist(intermediateFilename, 'file'))
          clear data;
          data = load(intermediateFilename);
          if(isfield(data,'summary'))
            surfFrameCount = (length(data.summary.surfResample) + 1) * frameSkip;
            while(surfFrameCount <= data.maxFrameCount)
              surfAsymp = asymptoteCoherenceReliability2(data.surfCoh, surfFrameCount);
              surfAsymp = rmfield(surfAsymp, 'percentiles');
              data.summary.surfResample(end+1) = surfAsymp;
              surfFrameCount = surfFrameCount + data.frameSkip;
            end
            icaFrameCount = (length(data.summary.icaResample) + 1) * frameSkip;
            while(icaFrameCount <= data.maxFrameCount)
              icaAsymp = asymptoteCoherenceReliability2(data.icaCoh, icaFrameCount);
              icaAsymp = rmfield(icaAsymp, 'percentiles');
              data.summary.icaResample(end+1) = icaAsymp;
              icaFrameCount = icaFrameCount + data.frameSkip;
            end
            data.summary.icaLabels = data.icaCoh.labels;
            data.summary.surfLabels = data.surfCoh.labels;
            summary = data.summary;
            path = fullfile(outputFolder, [file, '.mat']);
            save(path, 'summary', '-v7.3');
            delete(intermediateFilename);
            processThisFile = false;
          end
        else
          processThisFile = false;
        end
      end
    end
  end
  if(processThisFile)
    
    %load file and compute coherence
    if(veryShort)
      load('/home/data/EEG/processed/quickBdf.mat');
    else
      eeg = loadBdf(fileName);
    end
    
    summary = performReliabilityAnalysis2(eeg, outPath);
    %summary = performReliabilityAnalysis1(eeg, outPath);
    %summary = performReliabilityAnalysis(eeg, outPath);
    save(outPath, 'summary');
  end
end

if(doMaster)
  outFile = ['masterIca.mat'];
  outPath = fullfile(outputFolder, outFile);
  if(~exist(outPath, 'file'))
    placeHolder = sprintf('created on %s', char(datetime));
    save(outPath, 'placeHolder');
    allPaths = cell(length(fileNames), 1);
    for i = 1:length(fileNames)
      allPaths{i} = fullfile(inputFolder, fileNames(i).name);
    end
    icaCoh = deriveIcaCoherenceMatrix(allPaths, 8);
    save(outPath, icaCoh);
  end
end

if(false)
  clf;
  measureCount = size(summary.icaRhos, 2);
  ind = (measureCount - 31*5) : measureCount;
  plot(summary.sampleFrameDuration ./ (128 * 60), summary.icaRhos(:,ind));
  legend(summary.icaLabels(ind));
end

clear;