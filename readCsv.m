function [ t ] = readCsv( filename, tryToParseNumbers )
%READCSV reads a csv file in as a table
%   It is assumed that the first row corresponds to variable names.

newLine = sprintf('\r\n');
delim = ',';
escape = '"';

%read the file
rawData = fileread(filename);

%remove escaped newlines and delimeters
newlineInd = strfind(rawData, newLine);
commaInd = strfind(rawData, delim);
if(length(newlineInd) == 0)
    newLine = sprintf('\n');
    newlineInd = strfind(rawData, newLine);
end
escapedNewlineInd = false(size(newlineInd));
escapedCommaInd = false(size(commaInd));
inQuotes = false;
for i = 1:length(rawData)
    if(rawData(i) == '"')
        if(~inQuotes)
            inQuotes = true;
        else
            inQuotes = false;
        end
    elseif(rawData(i) == ',')
        if(inQuotes)
            escapedCommaInd(find(commaInd == i)) = true;
        end
    elseif(rawData(i) == sprintf('\n'))
        if(inQuotes)
            if(length(newLine) == 1)
                escapedNewlineInd(find(newlineInd == i)) = true;
            else
                escapedNewlineInd(find(newlineInd == i-1)) = true;
            end
        end
    end
end
commaInd(escapedCommaInd) = [];
newlineInd(escapedNewlineInd) = [];


%add beginning of file as fake newline
newlineInd = [1-length(newLine), newlineInd];
%add to end of file if the last line isn't blank
if(length(rawData) > max(newlineInd) + length(newLine) - 1)
    newlineInd = [newlineInd, length(rawData)+1];
end
t = table();
%iterate through each row
for i = 1:(length(newlineInd)-1)
    rowStart = newlineInd(i) + length(newLine) - 1;
    rowEnd = newlineInd(i+1);
    %take cells from within the row
    commas = commaInd((commaInd > rowStart) & (commaInd < rowEnd));
    commas = [rowStart, commas, rowEnd];
    for j = 1:(length(commas)-1)
        %fprintf('\ni=%d, j=%d', i, j);
        cellStart = commas(j)+1;
        cellEnd = commas(j+1)-1;
        if(cellEnd >= cellStart)
            item = rawData(cellStart:cellEnd);
            %remove unescaped quotes
            if(item(1) == '"' && item(end) == '"')
                item = item(2:end-1);
                %replace double quotes (escaped) with single (unescaped)
                item = strrep(item, '""', '"');
            end
        else
            item = '';
        end
        %first row contains column names
        if(i == 1)
            variableNames{j} = matlabSafeVariableName(item);
        %second row is used to determine if the column is numeric
        elseif(i == 2)
            number = str2double(item);
            if(length(number) > 0 && ~isnan(number))
                isNumeric(j) = 1;
                t{i-1, j} = number;
            else
                t{i-1, j} = {item};
            end
        %other rows simply supply data
        else
            if(isNumeric(j))
                number = str2double(item);
                t{i-1, j} = number;
            else
                t{i-1, j} = {item};
            end
        end
    end
    %perform initialization during early steps
    if(i == 1)
        isNumeric = false(size(variableNames));
    elseif(i == 2)
        t.Properties.VariableNames = variableNames;
        t{1:length(newlineInd) - 2,1} = t{1,1};
    end
end




end

