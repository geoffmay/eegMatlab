folder = 'C:\Users\Neuro\Documents\3dPrintBrains'

files = dir(folder);
files = files(cellfun(@length, strfind({files.name}, '.stl')) > 0);

fileNumber = 1;

fileId = fopen(fullfile(folder, files(fileNumber).name));
header = fscanf(fileId, '%c');
fclose(fileId);

newlines = strfind(header, sprintf('\n'));
newlines = [0 newlines length(header)+1];

vertices = strfind(header, 'vertex');

coords = NaN(length(vertices), 3);

% for i = 1:length(vertices)
vCounter= 1;
for i = 1:(length(newlines) - 1)
    fprintf('%d of %d\n', i, length(newlines));
    %     nextLine = newlines(min(find(newlines > vertices(i)))) - 1;
    %     line = header(vertices(i):nextLine);
    target = '    vertex';
    line = header((newlines(i)+1):(newlines(i+1)-1));
    if(length(line) > length(target) && strcmp(line(1:length(target)), target))
        items = strsplit(line, ' ');
        %     coords(i,1) = str2double(items{2});
        %     coords(i,2) = str2double(items{3});
        %     coords(i,3) = str2double(items{4});
        coords(vCounter,1) = fastParseDouble(items{3});
        coords(vCounter,2) = fastParseDouble(items{4});
        coords(vCounter,3) = fastParseDouble(items{5});
        vCounter = vCounter + 1;
    end
end


uniqueVertices = NaN(size(coords));
uVCounter = 1;
for i = 1:size(coords, 1)
    fprintf('%d of %d (%d uniques)\n', i, size(coords, 1), uVCounter)
    isUnique = true;
    for j = 1:(uVCounter - 1)
        isMatch = true;
        for k = 1:size(coords, 2)
            if(coords(i,k) ~= uniqueVertices(j,k))
                isMatch = false;
            end
        end
        if(isMatch)
            isUnique = false;
        end
    end
    if(isUnique)
        uniqueVertices(uVCounter,:) = coords(i,:);
        uVCounter = uVCounter + 1;
    end
end

end