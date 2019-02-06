folder = 'C:\Users\Neuro\Documents\MATLAB\processed\Robi\test';

mkdir(folder);
maxFiles = 1000;

for i = 1:maxFiles
    filePaths{i} = fullfile(folder, sprintf('%d.raw', i));
end

errorCount = 0;
concurrentFiles = 0;

for i = 1:maxFiles
    fileIds(i) = fopen(filePaths{i}, 'w');
    if(fileIds(i) ~= -1)
        fwrite(fileIds(i), i, 'double');
        concurrentFiles = concurrentFiles + 1;
    else
        errorCout = errorCount + 1;
    end
end

for i = 1:maxFiles
    if(fileIds(i) ~= -1)
        fclose(fileIds(i));
    end
end

errorCount
concurrentFiles
