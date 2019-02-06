function [tms] = findTmsBookends(data, voltageThreshold, frameThreshold)
%FINDTMSBOOKENDS Summary of this function goes here
%   Detailed explanation goes here
%identify start and finish of TMS
if(~exist('voltageThreshold', 'var'))
    voltageThreshold = 2500;
end
if(~exist('frameThreshold', 'var'))
    frameThreshold = 5050;
end


meanData = mean(abs(data),2);
supraThreshold = find(meanData > voltageThreshold);
redundantSuprathreshold = find(diff(supraThreshold) == 1) + 1;
tms.risingEdge = supraThreshold;
tms.fallingEdge = supraThreshold;
tms.risingEdge(redundantSuprathreshold) = [];
tms.fallingEdge(redundantSuprathreshold-1) = [];
longPause = find(diff(tms.risingEdge) > frameThreshold) + 1;
longPause2 = find(diff(tms.fallingEdge) > frameThreshold);
tms.tmsStart = tms.risingEdge([1; longPause]);
tms.tmsFinish = tms.fallingEdge([longPause2; length(tms.fallingEdge)]);
end

