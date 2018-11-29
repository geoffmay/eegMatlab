load('/home/data/EEG/data/Oregon/wahbehVariables.mat');

inputFolders = {'/home/data/EEG/data/Oregon/Auditory1', '/home/data/EEG/data/Oregon/Auditory2'};
outputFolder = '/home/data/EEG/processed/Oregon/phaseSlopeWindows';

if(~exist(outputFolder, 'file'))
  mkdir(outputFolder);
end

fileCounter = 1;
for i = 1:length(inputFolders)
  folder = inputFolders{i};
  files = dir(folder);
  files([files.isdir]) = [];
  for j = 1:length(files)
    inputFiles{fileCounter} = fullfile(folder, files(j).name);
    fileCounter = fileCounter + 1;
  end
end

for fileCounter = 1:length(inputFiles)
  inputFilename = inputFiles{fileCounter};
  [folder, file, ext] = fileparts(inputFilename);
  fprintf('\n%s: %s', char(datetime), file);
  eeg = loadBdf(inputFilename);
  keepChannels = 1:31;
  eeg.data = eeg.data(keepChannels,:);
  eegSlice = eeg;
  
  type = [eeg.event.type]';
  latency = [eeg.event.latency]';
  duration = diff([latency; eeg.pnts]);
  tab1 = table(latency, type, duration);
  tab2 = tabulate(type);
  tab2(tab2(:,2) == 0, :) = [];
  event.longform = tab1;
  event.summary = tab2;
  event.filename = inputFilename;
  events(fileCounter) = event;
  
  
  targetEventType = 44;
  preEvent = eeg.srate * -1;
  postEvent = eeg.srate * 2.5;
  targetIndices = find([eeg.event.type] == targetEventType);
  if(length(targetIndices) < 2)
    targetEventType = 1;
    targetIndices = find([eeg.event.type] == targetEventType);
  end
  clear topographies;
  for i = 1:length(targetIndices)
    latency = eeg.event(targetIndices(i)).latency;
    startIndex = latency + preEvent;
    endIndex = latency + postEvent;
    endIndex = min(endIndex, eeg.pnts);
    outputFile = sprintf('%s-%03d.mat', file, i);
    outputPath = fullfile(outputFolder, outputFile);
    if(~exist(outputPath, 'file'))
      fprintf('\n%s: %s', char(datetime), outputFile);
      stub = sprintf('started on %s', char(datetime));
      save(outputPath, 'stub');
      if(startIndex < 0)
        eegSlice.data = zeros(size(eeg.data,1), endIndex - startIndex + 1);
        destStart = -startIndex + 2;
        eegSlice.data(:, destStart:end) = eeg.data(:, 1:endIndex);
        for meanIndex = 1:size(eegSlice.data, 1)
          eegSlice.data(meanIndex, 1:(destStart-1)) = mean(eegSlice.data(meanIndex, destStart:end));
        end
        %debug
        toPlot = eegSlice.data;
        for meanIndex = 1:size(toPlot, 1)
          toPlot(meanIndex, :) = toPlot(meanIndex, :) - mean(toPlot(meanIndex, :));
        end
        plot(toPlot');
        pan xon;
        zoom xon;
        %end debug
      else
        eegSlice.data = eeg.data(:, startIndex:endIndex);
      end
      topography = phaseSlopeTimecourse1(eegSlice, [1 512]);
      if(startIndex < 0)
        topography.note = sprintf('filled sample 1-%d with zeros', startIndex);
      end
      topography.rootData = eegSlice.data;
      topography.details.preEventFrames = preEvent;
      topography.details.postEventFrames = postEvent;
      topography.details.targetEventType = targetEventType;
      
      save(outputPath, 'topography', '-v7.3');
    end
  end
end

save('/home/data/EEG/processed/Oregon/eventSummary.mat', 'events');


% inputFilename = '/home/data/EEG/data/Oregon/Auditory1/VM101.1.Tones.bdf';


%debug

chanIndex = 0;


close all;
x = (-1024:(1024 * 2.5)) ./ 1024;
chanIndex = chanIndex + 1;
y1 = topography.estimatedTimeLag(:, chanIndex);
y1 = y1 ./ nanstd(y1);
a = 1;
b0 = 10;
b = repmat(1/b0, [1, b0]);
y1f = filter(b, a, y1);
y2 = topography.rootData(chanIndex, :)';
y2 = y2 - mean(y2);
y2 = y2 ./ std(y2);
figure;
plot(x, [y1f y2]);
legend({'phaseSlope', 'rawData'});
title(sprintf('%s', topography.chanlocs(chanIndex).labels));
%end debug


% 
% 
% 
% startIndex = 174000;
% endIndex = 178000;
% 
% eeg.data = eeg.data(keepChannels, startIndex:endIndex);
% phaseSlopeTimecourse(eeg)
% 
% latency = ([eeg.event.latency] ./ eeg.srate)';
% type = [eeg.event.type]';
% duration = diff([latency; (eeg.pnts / eeg.srate)]);
% eventTab = table(latency, type, duration);

%eventTab(eventTab(:,2) == 0,:) = [];

