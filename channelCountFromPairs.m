function [ channelCount ] = channelCountFromPairs( channelPairCount )
%CHANNELCOUNTFROMPAIRS Calculates the number of channels from the number of
%channel pairs.
%   Since channelPairs = channels * (channels-1) / 2, this solves for
%   channels using the quadratic formula.

%pairs = (channels * (channels-1)) / 2
%0 = channels^2 - channels - 2*pairs;

%a = 1
%b = -1
%c = -2*pairs
%channels = (1 +/- sqrt(i - (4* -2* pairs))) / 2

channelCount = (1 + sqrt(1 + 8 * channelPairCount)) / 2;


end

