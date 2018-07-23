function result = countGoodMriTime(subjectNameFilter)

folder = '/home/data/subjects';
debug = false;

if(~exist('subjectNameFilter', 'var'))
  subjectNameFilter = 'TEMI';
end

subFolders = dir(folder);
subFolders(~[subFolders.isdir]) = [];
isTemi = cellfun(@length, strfind({subFolders.name}, subjectNameFilter)) > 0;
temiInd = find(isTemi);
fileCounts = zeros(length(temiInd), 3);
goodscanCounts = zeros(length(temiInd), 3);
totalscanCounts = zeros(length(temiInd), 3);

for i = 1:length(temiInd)
  if(debug)
    fprintf('\nreading subject %d of %d', i, length(temiInd));
  end
  filePath = fullfile(folder, subFolders(temiInd(i)).name);
  rawFolder = fullfile(filePath, 'raw');
  
  %load runs/sessions file
  %   sessionCatalogFile = fullfile(filePath, 'fc_processed', 'RSFC_runs_sessions.txt');
  fcProcessedFolder = fullfile(filePath, 'fc_processed');
  fcNames = dir(fcProcessedFolder);
  sessionCatalogFile = fcNames(find(cellfun(@length, strfind({fcNames.name}, 'runs_sessions'))>0));
  if(length(sessionCatalogFile)>1)
    sessionCatalogFile = sessionCatalogFile(1);
  end
  sessionCatalog = textread(fullfile(fcProcessedFolder, sessionCatalogFile.name));
  
  %load runs file
  %   tmaskFile = fullfile(filePath, 'fc_processed', 'RSFC_all_tmask.txt');
  tmaskFile = fcNames(find(cellfun(@length, strfind({fcNames.name}, 'tmask.txt'))>0));
  if(length(sessionCatalogFile)>1)
    tmaskFile = tmaskFile(1);
  end
  tmask = textread(fullfile(fcProcessedFolder, tmaskFile.name));
  
  %loop through all sessions
  sessions = unique(sessionCatalog(:,2));
  for j = 1:length(sessions)
    sessionInd = sessionCatalog(:,2) == sessions(j);
    totalCount = sum(sessionInd);
    goodCount = sum(tmask(sessionInd));
    if(sessions(j) < 100)
      goodscanCounts(i, 1) = goodscanCounts(i, 1) + goodCount;
      totalscanCounts(i, 1) = totalscanCounts(i, 1) + totalCount;
    elseif(sessions(j) < 200)
      goodscanCounts(i, 2) = goodscanCounts(i, 2) + goodCount;
      totalscanCounts(i, 2) = totalscanCounts(i, 2) + totalCount;
    else
      goodscanCounts(i, 3) = goodscanCounts(i, 3) + goodCount;
      totalscanCounts(i, 3) = totalscanCounts(i, 3) + totalCount;
    end
  end

  
  subSubFolders = dir(rawFolder);
  
  for j = 3:length(subSubFolders)
    scanName = subSubFolders(j).name;
    numberText = scanName(9:end);
    number = str2double(numberText);
    if(number < 100)
      fileCounts(i, 1) = fileCounts(i, 1) +1;
    elseif(number < 200)
      fileCounts(i, 2) = fileCounts(i, 2) +1;
    else
      fileCounts(i, 3) = fileCounts(i, 3) +1;
    end
  end
end

goodScanPercent = goodscanCounts ./ totalscanCounts .*100 ;
goodScanMinutes = goodscanCounts .* 3 ./ 60;
subjectId = {subFolders(temiInd).name}';

preMinutes = goodScanMinutes(:,1);
postMinutes = goodScanMinutes(:,2);
fuMinutes = goodScanMinutes(:,3);

pctPre = goodScanPercent(:,1);
pctPost = goodScanPercent(:,2);
pctFU = goodScanPercent(:,3);

result = table(subjectId, preMinutes, postMinutes, fuMinutes, pctPre, pctPost, pctFU);


