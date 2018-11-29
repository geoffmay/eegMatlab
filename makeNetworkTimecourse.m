subjectName = 'GeoffEEGTest';

%copy from external hard drive
copyfiles_frommedia_manualmount_BIDS(subjectName);

rootFolder = ['/home/data/subjects/sub-', subjectName, '/'];
infomapFolder = ['/home/data/subjects/sub-', subjectName, '/infomap/RSFC'];
infomap = ft_read_cifti_mod(fullfile(infomapFolder, 'rawassn_minsize100_regularized.dscalar.nii'));

sessionFolders = dir(rootFolder);
sessionFolders = sessionFolders(cellfun(@length, strfind({sessionFolders.name}, 'ses-')) > 0);
totalCounter = 1;
for sessionNumber = 1:length(sessionFolders)
  sessionFolder = fullfile(rootFolder, sessionFolders(sessionNumber).name);
  functionalFolder = fullfile(sessionFolder, 'func');
  if(exist(functionalFolder, 'file'))
    runFiles = dir(functionalFolder);
    runFiles = runFiles(cellfun(@length, strfind({runFiles.name}, '.json')) > 0);
    for runCounter = 1:length(runFiles)
      runFile = fullfile(functionalFolder, runFiles(runCounter).name);
      fileId = fopen(runFile);
      json = fscanf(fileId, '%c');
      fclose(fileId);
      jsonTarget = '"AcquisitionTime": "';
      startInd = strfind(json, jsonTarget) + length(jsonTarget);
      quotes = strfind(json, '"');
      endInd = min(quotes(quotes > startInd))-1;
      times{totalCounter} = json(startInd:endInd);
      items = strsplit(runFile, '_');
      isTask = find(cellfun(@length, strfind(items, 'task-')) > 0);
      items1 = strsplit(items{isTask}, '-');
      tasks{totalCounter} = items1{2};
      sessions{totalCounter} = sessionNumber;
      runs{totalCounter} = runCounter;
      runContentFile = strrep(runFile, '.json', '.nii.gz');
      runContent = load_untouch_nii_2D(runContentFile);
      volumeCounts{totalCounter} = size(runContent.img, 2);
      
      totalCounter = totalCounter + 1;
    end
  end
end

outputFilename = ['/home/data/Analysis/BOLD_EEG/', subjectName, '_metadata.csv'];
fileId = fopen(outputFilename, 'w');
fprintf(fileId, 'index,session,run,volumeCount,time,task');
for i = 1:length(times)
  fprintf(fileId, '\n%d,%d,%d,%d,%s,%s', i, sessions{i}, runs{i}, volumeCounts{i}, times{i}, tasks{i});
end
fclose(fileId);

% functionalFolder = fullfile(rootFolder, 'functional');
% functionalFiles = dir(functionalFolder);
% functionalFiles = functionalFiles(cellfun(@length, strfind({functionalFiles.name}, '.nii.gz')) > 0);
% 
