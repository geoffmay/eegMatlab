% inputFolder = '/home/data/EEG/processed/Oregon/artifactRemoved';
% outputFolder = '/home/data/EEG/processed/Oregon/artifactRemovedReliability';
inputFolder = '/home/data/EEG/processed/Oregon/artifactRemoved';
outputFolder = '/home/data/EEG/processed/Oregon/artifactRemovedReliability';

files = dir(inputFolder);

files([files.isdir]) = [];

for fileCounter = 1:length(files)
  inputFilePath = fullfile(inputFolder, files(fileCounter).name);
  outputFilePath = fullfile(outputFolder, files(fileCounter).name);
  
  if(~exist(outputFilePath, 'file'))
    load(inputFilePath);
    fprintf('\n%s', char(datetime)); tic;
    summary = performReliabilityAnalysis4(surfCoh, outputFilePath);
    fprintf('\n%s', char(datetime)); toc;
    save(outputFilePath, 'summary', '-v7.3');
  end
end

inputFolder = '/home/data/EEG/processed/Oregon/artifactPresent';
outputFolder = '/home/data/EEG/processed/Oregon/artifactPresentReliability';

files = dir(inputFolder);

files([files.isdir]) = [];

for fileCounter = 1:length(files)
  inputFilePath = fullfile(inputFolder, files(fileCounter).name);
  outputFilePath = fullfile(outputFolder, files(fileCounter).name);
  if(~exist(outputFilePath, 'file'))
    fprintf('\n%s', files(fileCounter).name);
    load(inputFilePath);
    [ surfCoh.matrix, surfCoh.labels ] = convertCoherenceStructToMatrix( surfCoh.coh, surfCoh.freqInfo, surfCoh.channels );
    fprintf('\n%s', char(datetime)); tic;
    summary = performReliabilityAnalysis4(surfCoh, outputFilePath);
    fprintf('\n%s', char(datetime)); toc;
    save(outputFilePath, 'summary', '-v7.3');
  end
end