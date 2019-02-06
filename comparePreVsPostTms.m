function [result] = comparePreVsPostTms(EEG)
%COMPAREPREVSresult.postTms Summary of this function goes here
%   Detailed explanation goes here

tms = findTmsBookends(EEG.data, 2500, EEG.sampleRate);
epochDuration = 1;
bufferDuration = 0.05;
epochLength = epochDuration * EEG.sampleRate;
bufferLength = bufferDuration * EEG.sampleRate;

%epoch x time x channel
result.preTms = NaN(length(tms.tmsStart), epochLength, EEG.chanCount);
result.postTms = NaN(length(tms.tmsStart), epochLength, EEG.chanCount);
preRange = NaN(length(tms.tmsStart), EEG.chanCount);
postRange = NaN(length(tms.tmsStart), EEG.chanCount);
rangeChange = NaN(length(tms.tmsStart), EEG.chanCount);
for i = 1:length(tms.tmsStart)
    startI = (tms.tmsStart(i) - epochLength - bufferLength + 1):(tms.tmsStart(i) - bufferLength);
    endI = (tms.tmsFinish(i)+bufferLength):(tms.tmsFinish(i) + epochLength + bufferLength - 1);
    result.preTms(i, :, :) = EEG.data(startI, :);
    result.postTms(i, :, :) = EEG.data(endI, :);
    for j = 1:EEG.chanCount
        preRange(i,j) = max(result.preTms(i,:,j)) - min(result.preTms(i,:,j));
        postRange(i,j) = max(result.postTms(i,:,j)) - min(result.postTms(i,:,j));
        rangeChange(i,j) = postRange(i,j) - preRange(i,j);
    end
end


%wavelet
for channelIndex = 1:EEG.chanCount
    preWav = NaN(size(result.preTms,1), 94, epochLength);
    postWav = NaN(size(result.preTms,1), 94, epochLength);
    for epochIndex = 1:size(result.preTms,1)
        fprintf('\nchannel %d, epoch %d', channelIndex, epochIndex);
        preData = squeeze(result.preTms(epochIndex, :, channelIndex));
        postData = squeeze(result.postTms(epochIndex, :, channelIndex));
        if(channelIndex == 1 && epochIndex == 1)            
            [waveletTransform, freqInfo, coneOfInfluence] = cwt(preData, EEG.sampleRate);
        end
        preWav(epochIndex, :, :) = cwt(preData);
        postWav(epochIndex, :, :) = cwt(postData);
    end
    %wavelet t test
    for freqCounter = 1:size(preWav, 2)
        fprintf('\nchannel %d freq %d of %d', channelIndex, freqCounter, size(preWav, 2));
        for timeCounter = 1:size(preWav, 3)
            x = preWav(:, freqCounter, timeCounter);
            y = postWav(:, freqCounter, timeCounter);
            [tTest.H,tTest.P,tTest.CI,tTest.STATS] = ttest(x,y);
            result.tTest.H(channelIndex, freqCounter, timeCounter) = tTest.H;
            result.tTest.P(channelIndex, freqCounter, timeCounter) = tTest.P;
            result.tTest.CiMin(channelIndex, freqCounter, timeCounter) = tTest.CI(1);
            result.tTest.CiMax(channelIndex, freqCounter, timeCounter) = tTest.CI(2);
            result.tTest.STATS(channelIndex, freqCounter, timeCounter) = tTest.STATS;
        end
    end
end

result.freqInfo = freqInfo;

% result.result.preTms = result.preTms;
% result.result.postTms = result.postTms;
% result.tTest = t;

%end

