function [ zScores ] = loadZScores( filename )
%LOADZSCORES Summary of this function goes here
%   Detailed explanation goes here

if(strcmp(filename((end-length('.zScorenotes.txt')+1):end), '.zScorenotes.txt'))
  headerFile = filename;
  dataFile = strrep(filename, '.zScorenotes.txt', '.zScore');
elseif(strcmp(filename(end-6:end), '.zScore'))
  dataFile = filename;
  headerFile = strrep(filename, '.zScore', '.zScorenotes.txt');
else
  error('unhandled file extension');
end

headerFileId = fopen(headerFile);
header = textread(headerFile,'%s');
%firstLineIndex = find(strcmp(header,'abs-Fp1-Delta1to4Hz,'));
firstLineIndex = find(strcmp(header,'time,'));
channelLabels = header(firstLineIndex:end);

endChunk = channelLabels(9:end);
startChunk = channelLabels(1:8);
channelLabels = [endChunk; startChunk];

fileInfo = dir(dataFile);
dataCount = fileInfo.bytes / 8;

columnCount = length(channelLabels);
rowCount = dataCount / columnCount;
rowCount = floor(rowCount);
dataCount = rowCount * columnCount;
fileId = fopen(dataFile);
%tic;
contents = fread(fileId, dataCount, 'single');
contents = single(contents);
%toc;
fclose(fileId);
data = reshape(contents, [columnCount, rowCount]);

zScores.timecourse = data;
zScores.labels = channelLabels;
% 
% sourceFreqs = {'-Delta1to4Hz,', '-Theta4to8Hz,', '-Alpha8to12Hz,', '-Beta12to25Hz,', '-HiBeta25to30Hz,'};
% destFreqs = {' delta', ' theta', ' alpha', ' beta', ' hibeta'};
% 
% goodFreq = zeros(size(data,1),1);
% for i = 1:length(sourceFreqs)
%     thisFreq = cellfun(@length, strfind(channelLabels, sourceFreqs{i}));
%     goodFreq = goodFreq | thisFreq;
% end
% 
% cohRow = cellfun(@length, strfind(channelLabels, 'coh'));
% absRow = cellfun(@length, strfind(channelLabels, 'abs'));
% relRow = cellfun(@length, strfind(channelLabels, 'rel'));
% 
% targetRows = (relRow | cohRow) & goodFreq;
% targetLabels = channelLabels(targetRows);
% 
% for i = 1:length(sourceFreqs)
%     targetLabels = strrep(targetLabels, sourceFreqs{i}, destFreqs{i});
% end
% targetLabels = strrep(targetLabels, 'coh-', '');
% for i = 1:length(targetLabels)
%     if(length(strfind(targetLabels{i}, 'rel-')) > 0)
%         targetLabels{i} = strrep(targetLabels{i}, 'rel-', '');
%         targetLabels{i} = [targetLabels{i}, ' rel'];
%     end
% end
% 
% plotData = mean(data(targetRows,:),2);
% keepData = range(data(targetRows,:),2) > 0;
% 
% values = plotData(keepData);
% labels = targetLabels(keepData);
% plotCoherencePca(plotData(keepData), targetLabels(keepData));
% 


end

