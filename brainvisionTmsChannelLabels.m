function labels = brainvisionTmsChannelLabels(headerFilename)

labels = cell(0);
headerFilename = 'C:\Users\Neuro\Downloads\TMS eeg\brainvision\patient001.vhdr';
fileId = fopen(headerFilename);
text = fscanf(fileId, '%c');
newlines = strfind(text, sprintf('\n'));

searchText = 'Resolution / Unit';
searchIndex = strfind(text, searchText);
firstLine = min(find(newlines > searchIndex));
searchText = 'S o f t w a r e  F i l t e r s';
searchIndex = strfind(text, searchText);
lastLine = min(find(newlines > searchIndex)) - 2;

for i = firstLine:lastLine
    line = text(newlines(i)+1:newlines(i+1)-2);
    items = strsplit(line, ' ');
    if(length(items) > 1)
        labels{end+1} = items{2};
    end
end

end