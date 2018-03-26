function [ channel1, channel2 ] = channelsFromPairIndex( pairIndex, maxPairs )
%CHANNELSFROMPAIRINDEX Summary of this function goes here
%   Detailed explanation goes here

channel1 = -1;
channel2 = -1;
maxChannel = channelCountFromPairs(maxPairs);

channelCounter = 1;
for i = 1:maxChannel
  for j = i+1:maxChannel
    if(channelCounter == pairIndex)
      channel1 = i;
      channel2 = j;
    end
    channelCounter = channelCounter + 1;
  end
end


end

