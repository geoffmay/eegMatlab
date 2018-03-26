function [ latency, height, area ] = maxErpSample( EEG, channelNumber )
%MAXERPSAMPLE Returns the 
%   Works best with a good band pass filter (e.g. 1-15 Hz).  Baselines
%   should also be removed.

allTrials = EEG.data(channelNumber,:,:);
voltSum = sum(allTrials, 3);
[maxVoltage, index] = max(voltSum);
latency = EEG.times(index);
[~,~,numberOfTrials] = size(allTrials);
height = maxVoltage / numberOfTrials;
area = 0;

end

