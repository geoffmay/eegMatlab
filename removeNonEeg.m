function [newEEG] = removeNonEeg(EEG)
% data = double.empty;
% channelLocations = EEG.chanlocs;
% channelLocations(1:length(channelLocations)) = [];
EEG = setStandardLocations(EEG);
newEEG = EEG;
newIndex = 1;
for i = 1:size(EEG.data,1)
    if(length(EEG.chanlocs(i).X) > 0)
%         data(size(data,1)+1,:,:) = EEG.data(i,:,:);
%         channelLocations(length(channelLocations)+1) = EEG.chanlocs(i);
        newIndex = newIndex + 1;
    else
        newEEG.nbchan = newEEG.nbchan-1;
        newEEG.data(newIndex, :) = [];
        newEEG.chanlocs(newIndex) = [];
    end
end


end