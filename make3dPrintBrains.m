inputFolder = '/home/data/subjects';
outputFolder = '/home/gmay/Documents/3dPrintBrains';

subFolders = dir(inputFolder);
for i = 3:length(subFolders)
  fprintf('\nfolder %d of %d (%s)', i, length(subFolders), char(datetime));
  subFolder = fullfile(inputFolder, subFolders(i).name);
  subFolder = [subFolder, '/fs_LR/NativeVol'];
  if(exist(subFolder, 'file'))
    files = dir(subFolder);
    isPial = cellfun(@length, strfind({files.name}, 'pial')) > 0;
    isSurf = cellfun(@length, strfind({files.name}, 'surf')) > 0;
    todo = find(isPial & isSurf);
    for j = 1:length(todo)
      inputFile = fullfile(subFolder, files(todo(j)).name);
      outputFile = fullfile(outputFolder, files(todo(j)).name);
      outputFile = strrep(outputFile, '.gii', '.stl');
      [a, b] = system(sprintf('mris_convert %s %s', inputFile, outputFile));
    end
  end
end