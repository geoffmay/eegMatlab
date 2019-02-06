

finishIncompleteFiles = 0;
makeDummy = 1; %creates a dummy file, which helps with running multiple instances
doIca = 0;



inputFolder = '/home/data/EEG/data/Oregon/PTSD';
outputFolder = '/home/data/EEG/processed/Oregon/coherence';

if(~exist(outputFolder))
  mkdir(outputFolder);
end
fileNames = dir(inputFolder);
fileNames([fileNames.isdir]) = [];


%loop through all files
for fileCounter = length(fileNames):-1:1
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
      end
    end
  end
  if(processThisFile)
    
    %load file and compute coherence
    eeg = loadBdf(fileName);
    eeg.data(32:end, :) = [];
    eeg.chanlocs(32:end) = [];
    eeg.nbchan = 31;
    %     end
    fprintf('\n%s: start ica', char(datetime));
    if(doIca)
      startTime = tic;
      [ica.weights,ica.sphere,ica.compvars,ica.bias,ica.signs,ica.lrates,ica.activations] = runica(eeg.data, 'maxSteps', 1024 * 16);
      summary.icaSecondsToComplete = toc(startTime);
      summary.ica = ica;
      save(outPath, 'summary');
    end
    
    [surfCoh, summary] = deriveRobiCoherenceMatrix(eeg);
    save(outPath, 'surfCoh', 'summary', '-v7.3');
  end
end
