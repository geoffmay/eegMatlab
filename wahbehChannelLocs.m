function [ labels, chanlocs ] = wahbehChannelLocs( input_args )
%WAHBEHCHANNELLOCS Summary of this function goes here
%   Detailed explanation goes here

labels = {'Fp1',   'AF3',   'F7',   'F3',   'FC1',   'FC5',   'T7',   'C3',   'CP1',   'CP5',   'P7',   'P3',   'Pz',   'PO3',   'O1',   'Oz',   'O2',   'PO4',   'P4',   'P8',   'CP6',...
    'CP2',   'C4',   'T8',   'FC6',   'FC2',   'F4',   'F8',   'AF4',   'Fp2',   'Fz',   'Cz',   'EXG1',   'EXG2',   'EXG3',   'EXG4',   'EXG5',   'EXG6',   'EXG7',   'EXG8',...
    'Ana1',   'Ana3',   'Status'};


if(nargout > 1)    
    
%     locs = load('C:\Users\Administrator\Documents\MATLAB\scripts\eeglab13_4_4b\chanlocs.mat');
    for i = 1:length(labels)
        EEG.chanlocs(i).labels = labels{i};
        
%         index = find(strcmp({locs.chanlocs.labels}, labels{i}));
%         if(length(index) > 0)
%             chanlocs(i) = locs.chanlocs(index);
%         end
    end
    EEG = setStandardLocations(EEG);
    chanlocs = EEG.chanlocs;
    
end

end

