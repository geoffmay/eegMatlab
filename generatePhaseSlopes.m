filenames = getRobiDataFiles();
outputFolder = '/home/data/EEG/processed/Robi/phaseSlopeReref/';

for i = 1:length(filenames)
  filename = filenames{i};
  [inputFolder,inputFilename,inputExtension] = fileparts(filename);
  outStart = strfind(inputFolder,'ROBI_');
  newOut = inputFolder(outStart:end);
  newOut = strrep(newOut, '/', '_');
  newOut = sprintf('%s_%s.mat', newOut, inputFilename);
  newOut = fullfile(outputFolder, newOut);

  outputFilename = sprintf('%s.mat',inputFilename);
  outputPath = fullfile(outputFolder,outputFilename);
  if(exist(outputPath,'file'))
    movefile(outputPath,newOut);
  elseif(~exist(newOut,'file'))
    now = char(datetime);
    fprintf('\n%s %s',now, filename);
    if(~exist('outputPath','file'))
      save(newOut, 'now');
      output = checkRereferencePhaseSlope(filename, false);
      save(newOut, 'output');
    end
    fprintf('(complete');
  end
end