function [ chanlocs ] = getChanlocs( labels )
%GETCHANLOCS Summary of this function goes here
%   Detailed explanation goes here

shortcut = '/home/data/EEG/processed/allLocs.mat'
if(exist(shortcut, 'file'))
  a = load(shortcut);
  chans = a.allLocs;
  clear a  
else
  
  path = fullfile(fileparts(which('eeglab')), 'plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
  chans = readlocs(path);
end

chanlocs = chans(1);
chanlocs(1) = [];

for i = 1:length(labels)
  for j = 1:length(chans);
    if(strcmp(lower(labels{i}), lower(chans(j).labels)))
      chanlocs(i) = chans(j);
    end
  end
end

end

