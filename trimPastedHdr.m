filename = 'pastedHdr.m';
fileId = fopen(filename);
text = fscanf(fileId, '%c');

fclose(fileId);
lineDelimeter = sprintf('\r\n');
newlines = strfind(text, lineDelimeter);
newlines = [(1-length(lineDelimeter)), newlines, (length(text) + 1)];
newText = '';
for i = 1:(length(newlines)-1)
    startI = newlines(i) + length(lineDelimeter);
    endI = newlines(i+1) - 1;
    line = text(startI:endI);
    if(length(line) > 2)
        line = line(3:end);
    else
        line = '';
    end
    newText = [newText, sprintf('\r\n'), line];
end