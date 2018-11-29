folders = {'/home/data/EEG/processed/Oregon/artifactPresentReliability', '/home/data/EEG/processed/Oregon/artifactRemovedReliability'};

inputPaths = cell(0);
for i = 1:length(folders)
  files = dir(folders{i});
  files([files.isdir]) = [];
  for j = 1:length(files)
    inputPaths{end+1} = fullfile(folders{i}, files(j).name);
  end
end

arraySize = 1;
allIterations = [];
for i = 1:length(inputPaths)
  fprintf('\nfile %d of %d', i, length(inputPaths));
  data = load(inputPaths{i});
  iterations = [data.summary.surfResample.iterations];
  allIterations(end+1:end+length(iterations)) = iterations;
end