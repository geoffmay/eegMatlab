function [ chanlocs ] = getChanlocs( labels )
%GETCHANLOCS Summary of this function goes here
%   Detailed explanation goes here

% path = fullfile(fileparts(which('eeglab')), 'plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
% chans = readlocs(path);
% 
% chanlocs = chans(1);
% chanlocs(1) = [];
% 
% for i = 1:length(labels)
%   for j = 1:length(chans);
%     if(strcmp(lower(labels{i}), lower(chans(j).labels)))
%       chanlocs(i) = chans(j);
%     end
%   end
% end

for chanCounter = 1:length(labels)
    EEG.chanlocs(chanCounter).labels = labels{chanCounter};;
end
EEG = setStandardLocations(EEG);
chanlocs = EEG.chanlocs;


end

