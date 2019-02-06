onlyEndpoints = 1;
makeDummy = 1;
files = getRobiDataFiles;

outputFolder = '/home/data/EEG/processed/Robi/coherence2';

if(onlyEndpoints)
  hasBase = cellfun(@length, strfind(files, 'baseline'))>0;
  hasExit = cellfun(@length, strfind(files, 'outcome'))>0;
  target = find(hasBase | hasExit);
else
  target = files;
end

for fileCounter = 1:length(target)
  filename = files{target(fileCounter)};
  slashes = strfind(filename, '/');
  sessionLabel = filename(slashes(end-2)+1:slashes(end)-1);
  sessionLabel =strrep(sessionLabel, '/', '-');
  outputPath = fullfile(outputFolder, sprintf('%s.mat', sessionLabel));
  if(~exist(outputPath, 'file'))
    if(makeDummy)
      placeHolder = sprintf('processing started on %s', char(datetime));
      save(outputPath, 'placeHolder');
    end
    thisFile.filename = filename;
    data = deriveRobiCoherenceMatrix(filename);
    if(onlyEndpoints)
      maxLength = 128 * 60 * 10;
    else
      if(length(strfind(filename, 'baseline')) > 0 || length(strfind(filename, 'outcome')) > 0)
        maxLength = 128 * 60 * 10;
      else
        maxLength = 128 * 60 * 30;
      end
    end
    if(size(data,1) > maxLength)
      data(maxLength+1:end,:) = [];
    end
    thisFile.means = mean(data.matrix,1);
    thisFile.stds = std(data.matrix,1);
    thisFile.skew = skewness(data.matrix,1);
    thisFile.kurtosis = kurtosis(data.matrix,1);
    save(outputPath, 'thisFile', '-v7.3');
  end
end