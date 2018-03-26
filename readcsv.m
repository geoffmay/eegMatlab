function [ t ] = readcsv( filename )
%READCSV Summary of this function goes here
%   Detailed explanation goes here

newLine = sprintf('\r\n');
delim = ',';
escape = '"';

%read the file
rawData = fileread(filename);
%split by newlines
rows = strsplit(rawData, newLine);
if(length(rows) == 1)
    rows = strsplit(rawData, sprintf('\n'));
    if(length(rows) == 1)
        rows = strsplit(rawData, sprintf('\r'));
    end
end
clear rawData;
%trim the newlines that strsplit didn't catch
for i = 1:length(rows)
    rows{i} = strrep(rows{i}, sprintf('\n'), '');
end
%split cells by delimeter
cells = cell(1, length(rows));
variableCount = 0;
for i = 1:length(rows)
    indices = strfind(rows{i}, delim);
    %ignore delimeters that are surrounded by escape characters
    escapes = strfind(rows{i}, escape);
    for j = 1:2:length(escapes)
        remove = find(escapes(j) < indices & indices < escapes(j+1));
        indices(remove) = [];
    end
    if(~isempty(indices))
        if(i ==1)
            variableCount = length(indices) + 1;
        elseif(length(indices) + 1 ~= variableCount)
            error(sprintf('invalid number of cells in row %d\r\n', i));
        end
        rowCells = cell(1, variableCount);
        text = rows{i};
        for j = 0:length(indices)
            if( j < length(indices))
                nextIndex = indices(j+1)-1;
            else
                nextIndex = length(rows{i});
            end
            if( j > 0)
                currentIndex = indices(j)+1;
            else
                currentIndex = 1;
            end
            if(currentIndex > nextIndex)
                rowCells{j+1} = '';
            else
                if(text(currentIndex) == escape)
                    if(text(nextIndex) == escape)
                        currentIndex = currentIndex + 1;
                        nextIndex = nextIndex - 1;
                    end
                end
                rowCells{j+1} = text(currentIndex:nextIndex);
            end
        end
        cells{i} = rowCells;
    else
        if(length(cells) >= i)
            cells(i) = [];
        end
    end
end
clear rows;
%initialize return value
t = table();
%this array will keep track of what's numeric to cut down on parsing time
isNumeric = zeros(1,variableCount);
%start at row 2, assuming first row is variable names
for j = 2:length(cells)
    row = cells{j};
    for i = 1:length(cells{j})
        number = 0;
        %figure out what's numeric on the first pass
        if(j == 2)
            number = str2double(row{i});
            if(length(number) > 0)
                 if(~isnan(number))
                    isNumeric(i) = 1;
                 elseif(strcmp(lower(row{i}),'nan')) 
                    isNumeric(i) = 1;
                 end
            end
        end
        %populate the table
        if(isNumeric(i))
            if(j > 2)
                number = str2double(row{i});
                if(length(number) == 0)
                    number = NaN;
                end
            end
            if(j ==5)
                dummy = 1;
            end
            t{j-1,i} = number;
        else
            t{j-1,i} = row(i);
        end
        %name the table variables on the first pass
        %(variables can't be named until the first element is added)
        if( j == 2)            
            topRow = cells{1};
            rowName = matlabSafeVariableName(topRow{i});
            newName = rowName;
            number = 0;
            while(sum(strcmp(t.Properties.VariableNames, newName)))
                number = number + 1;
                newName = strcat(rowName, num2str(number));
            end
            t.Properties.VariableNames{i} = newName;
        end
    end
    %allocate the rest of the table to prevent repeated allocations
    if(j==2)
        t{length(cells)-1,1} = t{1,1};
    end
        display(strcat('j:', num2str(j),' i:', num2str(i)));
end

end

function [ variableName ] = matlabSafeVariableName( variableName )
%MATLABSAFEVARIABLENAME replace unsafe characters
%   Replaces [#] [.] [$] [ ] respectively with 
%   [NUM] [POINT] [DOLLAR] [_]
            variableName = strrep(variableName, '#', 'NUM');
            variableName = strrep(variableName, '.', 'POINT');
            variableName = strrep(variableName, '$', 'DOLLAR');
            variableName = strrep(variableName, ' ', '_');
            variableName = strrep(variableName, '-', '_');
            variableName = strrep(variableName, '>', 'GreaterThan');
            variableName = strrep(variableName, '<', 'LessThan');
            
            if(length(variableName > 0))
            unsafestarts = {'_','0','1','2','3','4','5','6','7','8','9'};
            if(sum(strcmp(unsafestarts, variableName(1))))
                variableName = strcat('a', variableName);
            end
            else
                variableName = 'a';
            end
            if(length(variableName) > namelengthmax)
                variableName = variableName(1:namelengthmax);
            end
end