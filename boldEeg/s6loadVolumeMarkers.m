function [positions] = s6loadVolumeMarkers(markerName)
%LOADVOLUMEMARKERS Summary of this function goes here
%   Detailed explanation goes here
% fileId = fopen(markerName, 'r', 'native', 'UTF-8');
fileId = fopen(markerName);
text = fscanf(fileId, '%c');
fclose(fileId);
newlines = strfind(text, sprintf('\n'));
volumeMarker = '<Description>R128</Description>';
positionMarker = '<Position>';
positionEndMarker = '</Position>';
markerIndices = strfind(text, volumeMarker);
positionIndices = strfind(text, positionMarker);
positionEndIndices = strfind(text, positionEndMarker);
if(length(markerIndices) ~= length(positionIndices))
    keep = false(size(positionIndices));
    for i = 1:length(markerIndices)
        nextIndex = min(find(positionIndices > markerIndices(i)));
        keep(nextIndex) = true;
    end
    positionIndices = positionIndices(keep);
    positionEndIndices = positionEndIndices(keep);
%     error('file should contain only volume markers');
end

for i = 1:length(positionIndices)
    numStart = positionIndices(i) + length(positionMarker);
    numEnd = positionEndIndices(i)-1;
    numString = text(numStart:numEnd);
    positions(i) = str2double(numString);
end

end

