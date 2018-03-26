function [ antFiles1, antFiles2, tonesFiles1, tonesFiles2, ptsdFiles ] = sortFiles( varargin )
%SORTFILES Summary of this function goes here
%   Detailed explanation goes here



if(nargin == 0)    
    folder = '/Users/Geoff/Documents/EEG/VET MIND EEG files/'; %mac
    folder = '/home/data/EEG/data/';  %linux
else
    folder = varargin{1};
end

files= dir(folder);
filenames = {files.name};



%a = {files.name};

for k = length(filenames):-1:1
    if(nargin == 0)
        % remove folders
        if files(k).isdir
            files(k) = [ ];
            continue
        end
    end
    % remove folders starting with .
    fname = filenames{k};
    if fname(1) == '.'
        filenames(k) = [ ];
    end
    
end

b = strfind(filenames, 'ANT');

antFiles1 = [];
antFiles2 = [];
tonesFiles1 = [];
tonesFiles2 = [];
ptsdFiles = [];

searchStrings = [{'2.bdf'}, {'2_ANT.bdf'},{'2.ANT.bdf'}, {'_2_'}, {'.2.'}, {'T2'}];
for i = 1:length(filenames)
    if(strfind(filenames{i}, 'PM'))
        ptsdFiles = [ptsdFiles, filenames(i)];
    elseif(length(b{i})> 0) %ANT file
        found = false;
        for j = 1:length(searchStrings)
            if(length(strfind(filenames{i}, searchStrings{j})) > 0)
                found = true;
                antFiles2 = [antFiles2, filenames(i)];
                break;
            end
        end
        if(~found)
                antFiles1 = [antFiles1, filenames(i)];
        end
    else
                found = false;
        for j = 1:length(searchStrings)
            if(length(strfind(filenames{i}, searchStrings{j})) > 0)
                found = true;
                tonesFiles2 = [tonesFiles2, filenames(i)];
                break;
            end
        end
        if(~found)
                tonesFiles1 = [tonesFiles1, filenames(i)];
        end
    end
end


lists = [{antFiles1}, {antFiles2}, {tonesFiles1}, {tonesFiles2}];
folders = {'Flanker1', 'Flanker2', 'Auditory1', 'Auditory2'};

if(nargin == 0)
    for i = 1:length(lists)
        if(~exist(fullfile(folder, folders{i})))
            mkdir(fullfile(folder, folders{i}));
        end
        currentList = lists{i};
        for j = 1:length(currentList)
            source = fullfile(folder, currentList(j).name);
            destination = fullfile(folder, folders{i}, currentList(j).name);
            movefile(source, destination);
        end
    end
end

end

