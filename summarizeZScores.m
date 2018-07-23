function summarizeZScores(username, password)

inputFolder = '/External/MyBook/zScoresComplete/';
outputFolder = '/home/data/EEG/processed/Robi/zScoreSummaries/';

%list files on server
[code,fileBlock] = system(['curl -l -u ', username, ':', password, ' ftp://192.168.1.196' inputFolder]);
files = strsplit(fileBlock, sprintf('\n'));
hasZScoreText = cellfun(@length, strfind(files, '.zScore')) > 0;
isZScoreNotes = cellfun(@length, strfind(files, '.zScorenotes.txt')) > 0;
isZScoreData = hasZScoreText & ~isZScoreNotes;
zScoreDataIndexes = find(isZScoreData);

for fileCounter = 1:length(zScoreDataIndexes)
  %copy file and header to local machine
  filename = files{zScoreDataIndexes(fileCounter)};
  headername = strrep(filename, '.zScore', '.zScorenotes.txt');
  outputPath = [outputFolder filename '.mat'];
  if(~exist(outputPath, 'file'))
    
    fprintf('\n%s: processing file %d of %d (%s)', char(datetime), fileCounter, length(zScoreDataIndexes), filename);
    command = sprintf('curl -u %s:%s ftp://192.168.1.196%s%s -o %s%s', username, password, inputFolder, filename, outputFolder, filename);
    [code, message] = system(command);
    command = sprintf('curl -u %s:%s ftp://192.168.1.196%s%s -o %s%s', username, password, inputFolder, headername, outputFolder, headername);
    [code, message] = system(command);
    
    %process the data
    zScores = loadZScores([outputFolder filename]);
    zScoreSummary = summarizeMatrix(zScores.timecourse, 2, floor(size(zScores.timecourse,2)/8));
    zScoreSummary.labels = zScores.labels;
    zScoreSummary.filename = filename;
    save(outputPath, 'zScoreSummary', '-v7.3');
  end
  
  %delete local copy of input files
  delete([outputFolder, filename]);
  delete([outputFolder, headername]);
end
