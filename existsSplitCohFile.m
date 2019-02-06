function [ret] = existsSplitCohFile(edfFilename)
%TSSPLITCOHFILE Summary of this function goes here
%   Detailed explanation goes here
locations = fileLocations;
folder = fullfile(locations.eegFourier, edfFilename);
ret = exist(folder, 'file');

if(false)
    subdirs = dir(locations.eegFourier);
    subdirs(strcmp({subdirs.name}, '.')) = [];
    subdirs(strcmp({subdirs.name}, '..')) = [];
    for i = 1:length(subdirs)
        files(i) = length(dir(fullfile(locations.eegFourier, subdirs(i).name)));
        fprintf('%s: %d\n', subdirs(i).name, files(i));
    end
end

end

