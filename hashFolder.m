function [ total ] = hashFolder( ftpObj, folder, total )
%HASHFOLDER Summary of this function goes here
%   Detailed explanation goes here

skipAt = 1;

hashFoldername = '/Users/Geoff/Documents/MATLAB/filehash/temp';
opt.Input = 'File';
opt.Method = 'SHA-256';

hashFilename = '/Users/Geoff/Documents/MATLAB/filehash/hashValues.txt';
fileId = fopen(hashFilename, 'r');
hashContents = fscanf(fileId, '%c');
fclose(fileId);


contents = dir(ftpObj, folder);
if(skipAt)
    remove = false(size(contents));
    for i = 1:length(remove)
        filename = contents(i).name;
        if(filename(1) == '@')
            remove(i) = true;
        end
    end
    contents(remove) = [];
end
for i = 1:length(contents)
    if(~contents(i).isdir)
        source = fullfile(folder, contents(i).name);
        dest = fullfile(hashFoldername, contents(i).name);
        total = total + contents(i).bytes;
        fprintf('\n%d: %s', total, source);
        mget(ftpObj, source, dest);
        fullname = fullfile(hashFoldername, contents(i).name, folder, contents(i).name);
        hash = DataHash(fullname, opt);
        rmdir(fullfile(hashFoldername, contents(i).name), 's');
        logHash(hash, source, contents(i).bytes);
        
    else
        total = hashFolder(ftpObj, fullfile(folder, contents(i).name), total);
    end
end

end

function logHash(hash, filename, bytes)
hashFilename = '/Users/Geoff/Documents/MATLAB/filehash/hashValues.txt';
fileId = fopen(hashFilename, 'a');
fwrite(fileId, sprintf('%s|%s|%d\n', hash, filename, bytes));
fclose(fileId);
end

function result = isFileAlreadyHashed(hashData)

end
