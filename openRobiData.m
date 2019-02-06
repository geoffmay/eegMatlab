
%read data from csvs
if(~exist('data', 'var'))
    folder = 'C:\Users\Neuro\Documents\MATLAB\data\ROBI\neuropsych\CSV 2018-06-22';
    files = dir(folder);
    counter = 1;
    for i = 1:length(files)
        file = files(i).name;
        if(length(file) > 4)
            fprintf('\n%d of %d (%s)\n', i, length(files), file);
            if(strcmp(file(end-3:end), '.csv'))
                path = fullfile(folder, file);
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
                        if(contains(text, "ROBI"))
                            isMetaData(i) = false;
                            subjectIdColumnHeaders{i} = data(i).table.Properties.VariableNames{k};
                            fprintf("\n%s (%s)", subjectIdColumnHeaders{i}, data(i).filename);
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
cvltFilename = 'C:\Users\Neuro\Documents\MATLAB\data\ROBI\neuropsych\CVLT\ROBI_CVLT_expanded.csv';
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


outputFolder = 'C:\Users\Neuro\Documents\MATLAB\data\ROBI';
save(fullfile(outputFolder, 'neuropsychUnstructured.mat'), 'data', '-v7.3');
save(fullfile(outputFolder, 'neuropsych.mat'), 'masterTable', '-v7.3');

writeCsv('singleTable.csv', masterTable);


