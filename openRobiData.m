%function masterTable = openRobiData()


psychFolder = '/home/data/EEG/data/ROBI/psychTests/';
%read data from csvs
if(~exist('data', 'var'))
  %     folder = 'C:\Users\Neuro\Documents\MATLAB\data\ROBI\neuropsych\CSV 2018-06-22';
  
  accessExportFolder = fullfile(psychFolder, 'CSV 2018-06-22');
  files = dir(accessExportFolder);
  counter = 1;
  for i = 1:length(files)
    file = files(i).name;
    if(length(file) > 4)
      fprintf('\n%d of %d (%s)\n', i, length(files), file);
      if(strcmp(file(end-3:end), '.csv'))
        path = fullfile(accessExportFolder, file);
        datum.filename = file;
        datum.table = readcsv(path);
        data(counter) = datum;
        counter = counter + 1;
      end
    end
  end
end

%sort metaData
isMetaData = true(size(data));
subjectIdColumnHeaders = cell(size(data));
for i = 1:length(data)
  for j = 1:size(data(i).table, 1)
    if(isMetaData(i))
      for k = 1:size(data(i).table, 2)
        if(isMetaData(i))
          thisCell = data(i).table{j,k};
          if(iscell(thisCell))
            text = thisCell{1};
            if(length(strfind(text, 'ROBI')) > 0)
              isMetaData(i) = false;
              subjectIdColumnHeaders{i} = data(i).table.Properties.VariableNames{k};
              fprintf('\n%s (%s)', subjectIdColumnHeaders{i}, data(i).filename);
            end
          end
        end
      end
    end
  end
end

%sort tables into empty, raw, and meta data.
rawDataCounter = 1;
metaDataCounter = 1;
emptyTableCounter = 1;
for i = 1:length(isMetaData)
  if(strcmp(data(i).filename, 'ROBI Database Ohio General Quex.csv'))
    emptyTable(emptyTableCounter) = data(i);
    emptyTableCounter = emptyTableCounter + 1;
  elseif(strcmp(data(i).filename, 'Copy Of ROBI Database STROOP.csv'))
    emptyTable(emptyTableCounter) = data(i);
    emptyTableCounter = emptyTableCounter + 1;
  elseif(size(data(i).table, 1) == 0)
    emptyTable(emptyTableCounter) = data(i);
    emptyTableCounter = emptyTableCounter + 1;
  elseif(isMetaData(i))
    metaData(metaDataCounter) = data(i);
    metaDataCounter = metaDataCounter + 1;
  else
    rawData(rawDataCounter) = data(i);
    rawDataCounter = rawDataCounter + 1;
  end
end

%clean up subject names
subjectIds = cell(0);
for i = 1:length(rawData)
  ids = rawData(i).table{:,'Subj_ID'};
  for j = 1:length(ids)
    id = ids{j};
    newId = strrep(id, 'ROBI_', 'ROBI');
    targetLength = 7;
    if(length(newId) > targetLength)
      newId((targetLength+1):end) = [];
    end
    if(length(newId) < targetLength)
      numText = newId(5:end);
      number = str2double(numText);
      newId = [newId(1:4), sprintf('%03d', number)];
    end
    rawData(i).table{j,'Subj_ID'} = {newId};
  end
end

%find tables with multiple entries
isMultiEntry = false(size(rawData));
for i = 1:length(rawData)
  tab = rawData(i).table;
  sourceIds = rawData(i).table{:,'Subj_ID'};
  counts = tabulate(sourceIds);
  counts = [counts{:,2}];
  if(max(counts) > 1)
    isMultiEntry(i) = true;
    %debug
    fprintf('\n\n%s', rawData(i).filename);
    tabulate(sourceIds)
    %end debug
  end
end


%get unique ids
subjectIds = cell(0);
for i = 1:length(rawData)
  ids = rawData(i).table{:,'Subj_ID'};
  subjectIds = [subjectIds; ids];
end
Subj_ID = unique(subjectIds);

%stitch tables together
masterTable = table(Subj_ID);
for i = 1:length(rawData)
  tableName = rawData(i).filename;
  tableName = strrep(tableName, 'ROBI Database ', '');
  if(~isMultiEntry(i))
    tab = rawData(i).table;
  else
    %handle forms with multiple entries per subject differently
    sourceIds = rawData(i).table{:,'Subj_ID'};
    sourceTab = tabulate(sourceIds);
    sourceTab = sortrows(sourceTab, 1);
    counts = [sourceTab{:,2}];
    maxCount = max(counts);
    clear tab;
    tab = table();
    for j = 1:size(sourceTab, 1)
      id = sourceTab{j,1};
      tab{j,1} = {id};
      subTab = rawData(i).table(strcmp(sourceIds, sourceTab{j,1}), :);
      %convert to long form
      tabColumnCounter = 2;
      for k = 1:size(subTab, 1)
        for l = 1:size(subTab, 2)
          tab(j, tabColumnCounter) = subTab(k,l);
          tab.Properties.VariableNames{tabColumnCounter} = sprintf('%s_%d', subTab.Properties.VariableNames{l}, k);
          tabColumnCounter = tabColumnCounter + 1;
        end
      end
    end
    tab.Properties.VariableNames{1} = 'Subj_ID';
    
    tabColumnCounter = size(tab, 2) + 1;
    for j = 1:size(sourceTab, 1)
      id = sourceTab{j,1};
      subTab = rawData(i).table(strcmp(sourceIds, sourceTab{j,1}), :);
      %do table-specific processing
      if(strcmp(rawData(i).filename, 'Ohio_State_Step1_Baseline.csv'))
        numberOfInjuries = size(subTab, 1);
        tab{j, tabColumnCounter} = numberOfInjuries;
        blasts = subTab{:,'OS_Step1_Cause_Blast_Pre'} == 1;
        locNumber = subTab{:,'OS_Step1_LOC_Pre'};
        noLoc = locNumber == 1;
        mildLoc = locNumber == 2;
        modLoc = locNumber == 3;
        sevLoc = locNumber == 4;
        unkLoc = locNumber == 5;
        locNumber(locNumber == 5) = 0;
        severest = max(locNumber);
        severestBlast = max(locNumber(blasts));
        if(length(severestBlast) == 0)
          severestBlast = NaN;
        end
        severestNonBlast = max(locNumber(~blasts));
        if(length(severestNonBlast) == 0)
          severestNonBlast = NaN;
        end
        tab{j, tabColumnCounter+1} = severestBlast;
        tab{j, tabColumnCounter+2} = severestNonBlast;
        tab{j, tabColumnCounter+3} = severest;
        tab{j, tabColumnCounter+4} = sum(noLoc);
        tab{j, tabColumnCounter+5} = sum(mildLoc);
        tab{j, tabColumnCounter+6} = sum(modLoc);
        tab{j, tabColumnCounter+7} = sum(sevLoc);
      elseif(strcmp(rawData(i).filename, 'Ohio_State_Step3_Baseline.csv'))
        dummy = 1;
      end
    end
    if(strcmp(rawData(i).filename, 'Ohio_State_Step1_Baseline.csv'))
      tab.Properties.VariableNames{tabColumnCounter} = 'Number_Of_Injuries';
      tab.Properties.VariableNames{tabColumnCounter+1} = 'Severest_Blast';
      tab.Properties.VariableNames{tabColumnCounter+2} = 'Severest_NonBlast';
      tab.Properties.VariableNames{tabColumnCounter+3} = 'Severest';
      tab.Properties.VariableNames{tabColumnCounter+4} = 'NoLoc_Count';
      tab.Properties.VariableNames{tabColumnCounter+5} = 'Mild_Count';
      tab.Properties.VariableNames{tabColumnCounter+6} = 'Mod_Count';
      tab.Properties.VariableNames{tabColumnCounter+7} = 'Sev_Count';
    elseif(strcmp(rawData(i).filename, 'Ohio_State_Step3_Baseline.csv'))
    end
  end
  columnCounter = size(masterTable, 2) + 1;
  %add new table to master table
  for j = 1:size(tab, 2)
    varName = tab.Properties.VariableNames{j};
    if(~strcmp(varName, 'ID') && ~strcmp(varName, 'Subj_ID'))
      sampleCell = tab{1,j};
      if(isnumeric(sampleCell))
        nanColumn = table(NaN(size(masterTable,1),1));
        masterTable(:, columnCounter) = nanColumn;
      else
        cells = repmat({'null'}, size(masterTable,1),1);
        cellColumn = table(cells);
        masterTable(:, columnCounter) = cellColumn;
      end
      fullName = matlabSafeVariableName([tableName, varName]);
      sourceIds = tab{:,'Subj_ID'};
      for k = 1:size(masterTable,1)
        destId = Subj_ID{k};
        sourceRowIndex = find(strcmp(sourceIds, destId));
        if(length(sourceRowIndex) == 1)
          masterTable{k,columnCounter} = tab{sourceRowIndex, j};
        elseif(length(sourceRowIndex > 1))
          %todo: handle multiple hits
        end
      end
      masterTable.Properties.VariableNames{columnCounter} = fullName;
      columnCounter = columnCounter + 1;
    end
  end
  
end


%add CVLT data
cvltFilename = fullfile(psychFolder, 'CVLT', 'ROBI_CVLT_expanded.csv');
cvltData = readcsv(cvltFilename);
timepoints = {'base', 'exit', 'fu'};
for i = 1:length(timepoints)
  timepointColumnStart = size(masterTable, 2) + 1;
  timepointMatch = find(cellfun(@length, strfind(cvltData{:, 'inputFiles'}, timepoints{i})));
  for j = 1:length(timepointMatch)
    row = cvltData(timepointMatch(j), :);
    subjectId = row{1,1};
    subjectId = subjectId{1};
    subjectId = subjectId(1:7);
    destRowIndex = find(strcmp(masterTable{:,1}, subjectId));
    columnCounter = timepointColumnStart;
    for k = 2:size(row, 2)
      value = row{1, k};
      masterTable{destRowIndex, columnCounter} = value;
      if(j == 1)
        columnName = matlabSafeVariableName(sprintf('cvlt_%s_%s', timepoints{i}, row.Properties.VariableNames{k}));
        masterTable.Properties.VariableNames{columnCounter} = columnName;
      end
      columnCounter = columnCounter + 1;
    end
  end
end


%add EEG data
tic;
eegFolder = '/home/data/EEG/processed/Robi/summaries/';
eegFiles = dir(eegFolder);
eegFiles([eegFiles.isdir]) = [];
eegBaseline = cellfun(@length, strfind({eegFiles.name}, 'baseline eyes open'))>0;
eegBaseline = eegBaseline & ([eegFiles.bytes] > 4000000);
eegOutcome = cellfun(@length, strfind({eegFiles.name}, 'outcome eyes open'))>0;
eegOutcome = eegOutcome & ([eegFiles.bytes] > 3706514);

%put eeg measures (power, asymmetry, coherence) in the table
firstFile = true;
for eegFileCounter = 1:length(eegFiles)
  processThis = false;
  if(eegBaseline(eegFileCounter))
    timepoint = 'baseline';
    processThis = true;
  elseif(eegOutcome(eegFileCounter))
    timepoint = 'exit';
    processThis = true;
  end
  if(processThis)
    eegData = load(fullfile(eegFolder, eegFiles(eegFileCounter).name));
    if(firstFile)
      %allocate the memory in the table
      timepointCount = 3;
      baselineFirstColumn = size(masterTable, 2)+1;
      outcomeFirstColumn = size(masterTable, 2)+length(eegData.measureSummary.labels)+1;
      changeFirstColumn = size(masterTable, 2)+length(eegData.measureSummary.labels)*2+1;
      block = NaN(size(masterTable, 1), length(eegData.measureSummary.labels) * timepointCount);
      masterTable{:,end+1:end+size(block,2)} = block;
      %put titles in for each column
      for chanCounter = 1:length(eegData.measureSummary.labels)
        fprintf('\neeg column title %d of %d', chanCounter, length(eegData.measureSummary.labels));
        titleSuffix = eegData.measureSummary.labels{chanCounter};
        titleSuffix = strrep(titleSuffix, ' ', '_');
        titleSuffix = strrep(titleSuffix, '-', '_');
        baselineTitle = ['baseline_', titleSuffix];
        outcomeTitle = ['exit_', titleSuffix];
        changeTitle = ['change_', titleSuffix];
        
        %         baselineTitle = ['baseline_', eegData.measureSummary.labels{titleCounter}];
        %         baselineTitle = strrep(baselineTitle, ' ', '_');
        %         baselineTitle = strrep(baselineTitle, '-', '_');
        %baselineTitle = matlabSafeVariableName(baselineTitle);
        %         outcomeTitle = ['exit_', eegData.measureSummary.labels{titleCounter}];
        %         outcomeTitle = strrep(outcomeTitle, ' ', '_');
        %         outcomeTitle = strrep(outcomeTitle, '-', '_');
        %         changeTitle = ['change_', eegData.measureSummary.labels{titleCounter}];
        %         changeTitle = strrep(changeTitle, ' ', '_');
        %         changeTitle = strrep(changeTitle, '-', '_');
        %outcomeTitle = matlabSafeVariableName(outcomeTitle);
        masterTable.Properties.VariableNames{chanCounter+baselineFirstColumn-1} = baselineTitle;
        masterTable.Properties.VariableNames{chanCounter+outcomeFirstColumn-1} = outcomeTitle;
        masterTable.Properties.VariableNames{chanCounter+changeFirstColumn-1} = changeTitle;
      end
      firstFile = false;
    end
    subjectId = eegFiles(eegFileCounter).name(1:8);
    subjectId = strrep(subjectId, '_', '');
    dataRow = find(strcmp(masterTable{:,1}, subjectId));
    if(strcmp(timepoint, 'baseline'))
      firstColumn = baselineFirstColumn -1;
      otherColumn = outcomeFirstColumn -1;
    elseif(strcmp(timepoint, 'exit'))
      firstColumn = outcomeFirstColumn -1;
      otherColumn = baselineFirstColumn -1;
    end
    sample = masterTable{dataRow, firstColumn+1};
    otherSample = masterTable{dataRow, otherColumn+1};
    if(isnan(sample))
      if(~isnan(otherSample))
        nowComplete = true;
      else
        nowComplete = false;
      end
      for columnCounter = 1:length(eegData.measureSummary.labels)
        fprintf('\neeg file %d of %d, column %d of %d', eegFileCounter, length(eegFiles), columnCounter, length(eegData.measureSummary.labels));
        columnIndex = firstColumn + columnCounter;
        masterTable{dataRow, columnIndex} = eegData.measureSummary.meanValue(columnCounter);
        if(nowComplete)
          baseVal = masterTable{dataRow, baselineFirstColumn -1+columnCounter};
          exitVal = masterTable{dataRow, outcomeFirstColumn -1+columnCounter};
          changeVal = exitVal-baseVal;
          masterTable{dataRow, changeFirstColumn -1+columnCounter} = changeVal;
        end
      end
    end
  end
end
toc;

% %put phase slope data in the folder
% tic;
% eegFolder = '/home/data/EEG/processed/Robi/summaries/';
% eegFiles = dir(eegFolder);
% eegFiles([eegFiles.isdir]) = [];
% eegBaseline = cellfun(@length, strfind({eegFiles.name}, 'baseline eyes open'))>0;
% eegBaseline = eegBaseline & ([eegFiles.bytes] > 4000000);
% eegOutcome = cellfun(@length, strfind({eegFiles.name}, 'outcome eyes open'))>0;
% eegOutcome = eegOutcome & ([eegFiles.bytes] > 3706514);
% 
% firstFile = true;
% for eegFileCounter = 1:length(eegFiles)
%   processThis = false;
%   if(eegBaseline(eegFileCounter))
%     timepoint = 'baseline';
%     processThis = true;
%   elseif(eegOutcome(eegFileCounter))
%     timepoint = 'exit';
%     processThis = true;
%   end
%   if(processThis)
%     eegData = load(fullfile(eegFolder, eegFiles(eegFileCounter).name));
%     if(firstFile)
%       %allocate the memory in the table
%       timepointCount = 3;
%       baselineFirstColumn = size(masterTable, 2)+1;
%       outcomeFirstColumn = size(masterTable, 2)+length(eegData.measureSummary.labels)+1;
%       changeFirstColumn = size(masterTable, 2)+length(eegData.measureSummary.labels)*2+1;
%       block = NaN(size(masterTable, 1), length(eegData.measureSummary.labels) * timepointCount);
%       masterTable{:,end+1:end+size(block,2)} = block;
%       %put titles in for each column
%       for chanCounter = 1:length(eegData.measureSummary.labels)
%         fprintf('\neeg column title %d of %d', chanCounter, length(eegData.measureSummary.labels));
%         titleSuffix = eegData.measureSummary.labels{chanCounter};
%         titleSuffix = strrep(titleSuffix, ' ', '_');
%         titleSuffix = strrep(titleSuffix, '-', '_');
%         baselineTitle = ['baseline_', titleSuffix];
%         outcomeTitle = ['exit_', titleSuffix];
%         changeTitle = ['change_', titleSuffix];
%         
%         %         baselineTitle = ['baseline_', eegData.measureSummary.labels{titleCounter}];
%         %         baselineTitle = strrep(baselineTitle, ' ', '_');
%         %         baselineTitle = strrep(baselineTitle, '-', '_');
%         %baselineTitle = matlabSafeVariableName(baselineTitle);
%         %         outcomeTitle = ['exit_', eegData.measureSummary.labels{titleCounter}];
%         %         outcomeTitle = strrep(outcomeTitle, ' ', '_');
%         %         outcomeTitle = strrep(outcomeTitle, '-', '_');
%         %         changeTitle = ['change_', eegData.measureSummary.labels{titleCounter}];
%         %         changeTitle = strrep(changeTitle, ' ', '_');
%         %         changeTitle = strrep(changeTitle, '-', '_');
%         %outcomeTitle = matlabSafeVariableName(outcomeTitle);
%         masterTable.Properties.VariableNames{chanCounter+baselineFirstColumn-1} = baselineTitle;
%         masterTable.Properties.VariableNames{chanCounter+outcomeFirstColumn-1} = outcomeTitle;
%         masterTable.Properties.VariableNames{chanCounter+changeFirstColumn-1} = changeTitle;
%       end
%       firstFile = false;
%     end
%     subjectId = eegFiles(eegFileCounter).name(1:8);
%     subjectId = strrep(subjectId, '_', '');
%     dataRow = find(strcmp(masterTable{:,1}, subjectId));
%     if(strcmp(timepoint, 'baseline'))
%       firstColumn = baselineFirstColumn -1;
%       otherColumn = outcomeFirstColumn -1;
%     elseif(strcmp(timepoint, 'exit'))
%       firstColumn = outcomeFirstColumn -1;
%       otherColumn = baselineFirstColumn -1;
%     end
%     sample = masterTable{dataRow, firstColumn+1};
%     otherSample = masterTable{dataRow, otherColumn+1};
%     if(isnan(sample))
%       if(~isnan(otherSample))
%         nowComplete = true;
%       else
%         nowComplete = false;
%       end
%       for columnCounter = 1:length(eegData.measureSummary.labels)
%         fprintf('\neeg file %d of %d, column %d of %d', eegFileCounter, length(eegFiles), columnCounter, length(eegData.measureSummary.labels));
%         columnIndex = firstColumn + columnCounter;
%         masterTable{dataRow, columnIndex} = eegData.measureSummary.meanValue(columnCounter);
%         if(nowComplete)
%           baseVal = masterTable{dataRow, baselineFirstColumn -1+columnCounter};
%           exitVal = masterTable{dataRow, outcomeFirstColumn -1+columnCounter};
%           changeVal = exitVal-baseVal;
%           masterTable{dataRow, changeFirstColumn -1+columnCounter} = changeVal;
%         end
%       end
%     end
%   end
% end
% toc;

%put phase slope measures in the table
% [chanlabels, chanlocs] = antChannelLocs;
% chanlocs = chanlocs(1:33);
firstFile = true;
zScoreFolder = '/home/data/EEG/processed/Robi/zScoreSummaries';
zScoreFiles = dir(zScoreFolder);
zScoreBaseline = cellfun(@length, strfind({zScoreFiles.name}, 'baseline_eyes_open'))>0;
zScoreOutcome = cellfun(@length, strfind({zScoreFiles.name}, 'outcome_eyes_open'))>0;

firstSession = false(size(zScoreFiles));
lastSession = false(size(zScoreFiles));
for i = 1:100
  target = sprintf('ROBI_%03d-tx_', i);
  match = find(cellfun(@length, strfind({zScoreFiles.name}, target)));
  if(length(match) > 0)
    maxInd = -1; 
    maxValue = realmin;
    minInd = -1;
    minValue = realmax;
    for j = 1:length(match)
      ind = match(j);
      name = zScoreFiles(ind).name;
      numberText = name((length(target)+1):end);
      dashInd = strfind(numberText, '-');
      numberText(dashInd:end) = [];
      number = str2num(numberText);
      if(number > maxValue)
        maxInd = ind;
        maxValue = number;
      end
      if(number < minValue)
        minInd = ind;
        minValue = number;
      end
    end
    if(minInd > -1)
      firstSession(minInd) = true;
    end
    if(maxInd > -1)
      lastSession(maxInd) = true;
    end
  end
end

%todo: use first and last sessions to populate table.

for zScoreFileCounter = 1:length(zScoreFiles)
  processThis = false;
  if(zScoreBaseline(zScoreFileCounter))
    timepoint = 'baseline';
    processThis = true;
  elseif(zScoreOutcome(zScoreFileCounter))
    timepoint = 'exit';
    processThis = true;
  end
  if(processThis)
    zScoreData = load(fullfile(zScoreFolder, zScoreFiles(zScoreFileCounter).name));
    if(firstFile)
      %allocate the memory in the table
      originalTableSize = size(masterTable, 2);
      chanCount = length(zScoreData.zScoreSummary.meanValue);
      block = NaN(size(masterTable, 1), chanCount * 3);
      baselineFirstColumn = originalTableSize + 1;
      outcomeFirstColumn = originalTableSize + 1 + chanCount;
      changeFirstColumn = originalTableSize + 1 + chanCount * 2;
      masterTable{:,end+1:end+size(block,2)} = block;
      %set titles
      for chanCounter = 1:length(zScoreData.zScoreSummary.labels)
        fprintf('\nzScore column title %d of %d', chanCounter, length(zScoreData.zScoreSummary.labels));
        %         titleSuffix = zScoreData.measureSummary.labels{chanCounter};
        titleSuffix = zScoreData.zScoreSummary.labels{chanCounter};
        titleSuffix = strrep(titleSuffix, ' ', '_');
        titleSuffix = strrep(titleSuffix, '-', '_');
        baselineTitle = ['zScore_baseline_', titleSuffix];
        outcomeTitle = ['zScore_exit_', titleSuffix];
        changeTitle = ['zScore_change_', titleSuffix];
        baselineTitle = matlabSafeVariableName(baselineTitle);
        outcomeTitle = matlabSafeVariableName(outcomeTitle);
        changeTitle = matlabSafeVariableName(changeTitle);
        
        masterTable.Properties.VariableNames{chanCounter+baselineFirstColumn-1} = baselineTitle;
        masterTable.Properties.VariableNames{chanCounter+outcomeFirstColumn-1} = outcomeTitle;
        masterTable.Properties.VariableNames{chanCounter+changeFirstColumn-1} = changeTitle;
      end
      firstFile = false;
    end
    %find destination row
    subjectId = zScoreFiles(zScoreFileCounter).name(1:8);
    subjectId = strrep(subjectId, '_', '');
    dataRow = find(strcmp(masterTable{:,1}, subjectId));
    %find destination column
    clear firstColumn otherColumn;
    if(strcmp(timepoint, 'baseline'))
      firstColumn = baselineFirstColumn -1;
      otherColumn = outcomeFirstColumn -1;
    elseif(strcmp(timepoint, 'exit'))
      firstColumn = outcomeFirstColumn -1;
      otherColumn = baselineFirstColumn -1;
    end
    sample = masterTable{dataRow, firstColumn+1};
    otherSample = masterTable{dataRow, otherColumn+1};
    if(isnan(sample))
      if(~isnan(otherSample))
        nowComplete = true;
      else
        nowComplete = false;
      end
      %copy data to table
      for columnCounter = 1:length(zScoreData.zScoreSummary.labels)
        fprintf('\nzScore file %d of %d, column %d of %d', zScoreFileCounter, length(zScoreFiles), columnCounter, length(zScoreData.zScoreSummary.labels));
        columnIndex = firstColumn + columnCounter;
        masterTable{dataRow, columnIndex} = zScoreData.zScoreSummary.meanValue(columnCounter);
        if(nowComplete)
          baseVal = masterTable{dataRow, baselineFirstColumn -1+columnCounter};
          exitVal = masterTable{dataRow, outcomeFirstColumn -1+columnCounter};
          changeVal = exitVal-baseVal;
          masterTable{dataRow, changeFirstColumn -1+columnCounter} = changeVal;
        end
      end
    end
  end
end




%this makes doing brute force correlations easier
masterTable = makeTableNumeric(masterTable);

%get some mri datas
% mriFile = '/home/data/Analysis/ROBI/DMN_deltas.txt';
% mriFileId = fopen(mriFile);
% mriText = fscanf(mriFileId, '%c');
% fclose(mriFileId);
% mriNewlines = strfind(mriText, sprintf('\n'));
% mriNewlines = [0 mriNewlines length(mriText)+1];
% mriColumnIndex = size(masterTable, 2) + 1;
% mriColumn = NaN(size(masterTable, 1), 1);
% for i = 1:(length(mriNewlines)-1)
%   line = mriText(mriNewlines(i)+1:mriNewlines(i+1)-1);
%   if(length(line) > 0)
%     subjectId = line(1:7);
%     rowIndex = find(strcmp(masterTable{:,1}, subjectId));
%     numberString = line(10:end);
%     number = str2double(numberString);
%     mriColumn(rowIndex) = number;
%   end
% end
% masterTable{:,mriColumnIndex} = mriColumn;
% masterTable.Properties.VariableNames{mriColumnIndex} = 'mri_dmn_delta';
mriFile = '/home/data/Analysis/ROBI/within_dmn.mat';
mriTab = load(mriFile);
mriTab = mriTab.tab;

mriColumnIndex = size(masterTable, 2) + 1;
mriAppend = NaN(size(masterTable, 1), size(mriTab, 2)-1);
for i = 1:(size(mriTab,1))
  if(length(line) > 0)
    subjectId = mriTab{i,1};
    rowIndex = find(strcmp(masterTable{:,1}, subjectId));
    mriAppend(rowIndex, :) = mriTab{i,2:end};
  end
end
columnRange = mriColumnIndex:mriColumnIndex+size(mriAppend,2)-1;
masterTable{:,columnRange} = mriAppend;
columnNames = {mriTab.Properties.VariableNames{2:end}};
masterTable.Properties.VariableNames(columnRange) = columnNames;


%debug
find(cellfun(@length, strfind(masterTable.Properties.VariableNames, 'cvlt')))
%end debug

%debug
% stroopInd = find(cellfun(@length, strfind({rawData.filename}, 'STROOP')));
% stroopTables = {rawData.filename};
% stroopTables = stroopTables(stroopInd);
% origStroop = rawData(find(strcmp({rawData.filename}, 'ROBI Database STROOP.csv')));
% copyStroop = rawData(find(strcmp({rawData.filename}, 'Copy Of ROBI Database STROOP.csv')));

% masterTable.Properties.VariableNames(find(cellfun(@length, strfind(lower(masterTable.Properties.VariableNames), 'age'))))'
% masterTable.Properties.VariableNames(find(cellfun(@length, strfind(lower(masterTable.Properties.VariableNames), 'dob'))))'
% masterTable.Properties.VariableNames(find(cellfun(@length, strfind(lower(masterTable.Properties.VariableNames), 'birth'))))'

%
% falseDuplicate = 'ROBI Database Tx Pt Hx Ohio Vasterling.csv';
% falseInd = find(cellfun(@length, strfind({rawData.filename}, falseDuplicate)));
% falseTable = rawData(falseInd);
% row1 = falseTable.table(2,:);
% row2 = falseTable.table(3,:);


%convert to long form
longformTable = table();
baselineColumns = cellfun(@length, strfind(masterTable.Properties.VariableNames, 'baseline_'))>0;
exitColumns = cellfun(@length, strfind(masterTable.Properties.VariableNames, 'exit_'))>0;
changeColumns = cellfun(@length, strfind(masterTable.Properties.VariableNames, 'change_'))>0;
baselineIndex = find(baselineColumns);
longValues = NaN(2*size(masterTable,1), length(baselineIndex));
for i = 1:length(baselineIndex)
  fprintf('\n%d of %d', i, length(baselineIndex));
  columnIndex = baselineIndex(i);
  columnTitle = masterTable.Properties.VariableNames{baselineIndex};
  exitTitle = strrep(columnTitle, 'baseline_', 'exit_');
  exitIndex = find(strcmp(masterTable.Properties.VariableNames, exitTitle));
  longValues(:, i) = [masterTable{:,columnIndex}; masterTable{:,exitIndex}];
end
longformTable{:,1} = [masterTable{:,1}; masterTable{:,1}];
longformTable{:,2:size(longValues,2)+1} = longValues;
longformTable.Properties.VariableNames(2:end) = strrep({masterTable.Properties.VariableNames{baselineIndex}}, 'baseline_', '');
longformTable.Properties.VariableNames(1) = masterTable.Properties.VariableNames(1);





outputFolder = psychFolder;
save(fullfile(outputFolder, 'neuropsychUnstructured.mat'), 'data', '-v7.3');
save(fullfile(outputFolder, 'neuropsych.mat'), 'masterTable', '-v7.3');
save(fullfile(outputFolder, 'neuropsychPrePostCombined.mat'), 'longformTable', '-v7.3');



%writeCsv(fullfile(outputFolder, 'neuropsych.csv'), masterTable);




