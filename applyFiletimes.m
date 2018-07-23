import java.io.File java.text.SimpleDateFormat


root = fileparts(which('preserveFiletimes'));

fileTimeFile = which('fileTimes');
fileId = fopen(fileTimeFile);
text = fscanf(fileId, '%c');
fclose(fileId);
newlines = strfind(text, sprintf('\n'));
newlines = [0 newlines length(text) + 1];

for i = 2:length(newlines)
    %parse date
    line = text(newlines(i-1)+1:newlines(i)-1);
    spaces = strfind(line, ' ');
    if(length(spaces) > 0)
        numberText = line(1:(spaces(1)-1));
        filePath = line(spaces(1)+1:end);
        number = str2double(numberText);
        filedate = datetime(number, 'convertfrom', 'datenum');
        
        %set file time (using java import)
        f = File(fullfile(root, filePath));
        sdf = SimpleDateFormat('dd-MMM-yyyy HH:mm:ss');
        newDate = sdf.parse(datestr(filedate));
        f.setLastModified(newDate.getTime);
    end    
end