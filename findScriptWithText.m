function [output] = findScriptWithText(target, root)
%FINDSCRIPTWITHTEXT Summary of this function goes here
%   Detailed explanation goes here
files = cell(0);
if(~exist('root', 'var'))
    root = 'C:\';
end
output = recursiveSearch(target, root, files);




function files = recursiveSearch(target, folder, files)

contents = dir(folder);
contents(strcmp({contents.name}, '.')) = [];
contents(strcmp({contents.name}, '..')) = [];

subFiles = contents(~[contents.isdir]);
for i = 1:length(subFiles)
    if(length(subFiles(i).name) > 3)
        if(strcmp(subFiles(i).name(end-1:end), '.m'))
            thisFile = fullfile(folder, subFiles(i).name);
            fileId = fopen(thisFile);
            text = fscanf(fileId, '%c');
            fclose(fileId);
            if(strfind(text, target))
                files{end+1} = thisFile;
                fprintf('%s\n', thisFile);
            end            
        end
    end
end
subFolders = contents([contents.isdir]);
for i = 1:length(subFolders)
    files = recursiveSearch(target, fullfile(folder, subFolders(i).name), files);
end
