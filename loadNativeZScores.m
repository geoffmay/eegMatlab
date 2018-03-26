folder = '/home/data/EEG/processed/Robi/zScoresComplete';
allFiles = dir(folder);
headerFile = 'ROBI_003-outcome_eyes_open-630230539149591228.zScorenotes.txt';
dataFile = 'ROBI_003-outcome_eyes_open-630230539149591228.zScore';

% headerFile = 'ROBI_003-baseline_eyes_open-630158995692243270.zScorenotes.txt';
% dataFile = 'ROBI_003-baseline_eyes_open-630158995692243270.zScore';

header = textread(fullfile(folder,headerFile),'%s');
%firstLineIndex = find(strcmp(header,'abs-Fp1-Delta1to4Hz,'));
firstLineIndex = find(strcmp(header,'time,'));
channelLabels = header(firstLineIndex:end);

endChunk = channelLabels(9:end);
startChunk = channelLabels(1:8);
channelLabels = [endChunk; startChunk];

fileInfo = dir(fullfile(folder,dataFile));
dataCount = fileInfo.bytes / 8;

columnCount = length(channelLabels);
rowCount = dataCount / columnCount;
rowCount = floor(rowCount);
dataCount = rowCount * columnCount;
fileId = fopen(fullfile(folder, dataFile));
tic;
contents = fread(fileId, dataCount, 'single');
contents = single(contents);
toc;
fclose(fileId);
data = reshape(contents, [columnCount, rowCount]);

sourceFreqs = {'-Delta1to4Hz,', '-Theta4to8Hz,', '-Alpha8to12Hz,', '-Beta12to25Hz,', '-HiBeta25to30Hz,'};
destFreqs = {' delta', ' theta', ' alpha', ' beta', ' hibeta'};

goodFreq = zeros(size(data,1),1);
for i = 1:length(sourceFreqs)
    thisFreq = cellfun(@length, strfind(channelLabels, sourceFreqs{i}));
    goodFreq = goodFreq | thisFreq;
end

cohRow = cellfun(@length, strfind(channelLabels, 'coh'));
absRow = cellfun(@length, strfind(channelLabels, 'abs'));
relRow = cellfun(@length, strfind(channelLabels, 'rel'));

targetRows = (relRow | cohRow) & goodFreq;
targetLabels = channelLabels(targetRows);

for i = 1:length(sourceFreqs)
    targetLabels = strrep(targetLabels, sourceFreqs{i}, destFreqs{i});
end
targetLabels = strrep(targetLabels, 'coh-', '');
for i = 1:length(targetLabels)
    if(length(strfind(targetLabels{i}, 'rel-')) > 0)
        targetLabels{i} = strrep(targetLabels{i}, 'rel-', '');
        targetLabels{i} = [targetLabels{i}, ' rel'];
    end
end

plotData = mean(data(targetRows,:),2);
keepData = range(data(targetRows,:),2) > 0;

values = plotData(keepData);
labels = targetLabels(keepData);
plotCoherencePca(plotData(keepData), targetLabels(keepData));

% %debug
% rowCount = 100;
% hex = cell(size(data,1), rowCount);
% for i = 1:size(data,1)
%     fprintf('\n%d',i);
%     for j = 1:rowCount
%         hex{i,j} = num2hex(data(i,j));
%     end
% end
% 
% isbig = data(:,rowCount) > 10;
% 
% hex2dec(hex(:,1))

% plot(dec);
% plot(contents(1:15000));
% pan xon;
% zoom xon;
% %end debug