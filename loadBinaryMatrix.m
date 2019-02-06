function data = loadBinaryMatrix(filename, numberOfColumns)

fileInfo = dir(filename);
fileLength = fileInfo.bytes / 8;
sampleCount = fileLength / numberOfColumns;
if(floor(sampleCount) ~= sampleCount)
  sampleCount = floor(sampleCount);
  fileLength = sampleCount * numberOfColumns;
end
% [folder, file, ext] = fileparts(filename);
% fileTimestamp = str2num(file);
% if(fileTimestamp > watershed)
%   multiplier = 1;
% end
fileId = fopen(filename);
contents = fread(fileId, fileLength, 'double');
fclose(fileId);

%truncate if necessary
% sampleCount = fileLength / numberOfColumns;
% if(sampleCount ~= floor(sampleCount))
%   sampleCount = floor(sampleCount);
%   fileLength = sampleCount * numberOfColumns;
%   contents = contents(1:fileLength);
% end
% labels = antChannelLocs;

data = reshape(contents, numberOfColumns, fileLength / numberOfColumns)';
% clear contents;

% if(rereference)
%   cpzIndex = find(strcmp(labels,'CPz'));
%   m1Index = find(strcmp(labels,'M1'));
%   m2Index = find(strcmp(labels,'M2'));
%   
%   for i=1:size(data,1)
%     sample = data(i,:);
%     sample = sample .* multiplier;
%     linkedMastoids = (sample(m1Index) + sample(m2Index)) / 2;
%     %     avg = mean(sample(1:33));
%     newSample = sample - linkedMastoids;
%     %     avgSample = sample - avg;
%     data(i,1:33) = newSample(1:33);
%   end
% end


end