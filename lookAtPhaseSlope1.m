folder = '/home/data/EEG/processed/Oregon/phaseSlopeWindows';
%folder = '/home/data/EEG/processed/Oregon/phaseSlope';

doPlot = false;
doNeuropsychCorrelations = false;
files = dir(folder);
singleSubject = true;

files(cellfun(@length, strfind({files.name}, 'ANT')) > 0) = [];

clear allEvents meanEvents;
psyData = load('/home/data/EEG/data/Oregon/wahbehVariables.mat');
originalColumnCount = size(psyData.vetmindData, 2);

if(singleSubject)
  subjects = [101];
else
  subjectNumbers = NaN(size(files));
  for i = 1:length(files)
    filename = files(i).name;
    if(length(filename) > 4)
      if(strcmp(filename(1:2), 'VM'))
        numberText = filename(3:5);
        subjectNumbers(i) = str2double(numberText);
      end
    end
  end
  subjects = unique(subjectNumbers);
  subjects(isnan(subjects)) = [];
end

for subjectCounter = 1:length(subjects)
  fprintf('\nsubject %d of %d', subjectCounter, length(subjects));
  %target = 'VM101';
  target = sprintf('VM%d', subjects(subjectCounter));
  hasTarget = cellfun(@length, strfind({files.name}, target)) > 0;
  targetFiles = files(hasTarget);
  clear eventData;
  for fileCounter = 1:length(targetFiles)
    data = load(fullfile(folder, targetFiles(fileCounter).name));
    if(isfield(data, 'topography'))
      if((subjectCounter == 1) && (fileCounter == 1))
        block = NaN(size(psyData.vetmindData, 1), size(data.topography.estimatedTimeLag, 2));
        psyData.vetmindData{:, (end+1):(end+size(block,2))} = block;
      end
      
      if(length(size(data.topography.estimatedTimeLag)) == 3)
        eventData(:,:,fileCounter) = data.topography.estimatedTimeLag(:,:,1);
      else
        eventData(:,:,fileCounter) = data.topography.estimatedTimeLag(:,:);
      end
    end
  end
  meanEvent = mean(eventData, 3);
  x = (((1:size(meanEvent,1)) - 1024) ./ 1024)';
  filteredEvent = eventData;
  filterDistance = 21;
  for i = 1:size(eventData,1)
    i1 = max(1, i-filterDistance);
    i2 = min(size(meanEvent,1), i+filterDistance);
    filteredEvent(i, :, :) = mean(eventData(i1:i2, :, :), 1);
  end
  endIndex = size(meanEvent, 2);
  %   endIndex = 1;
  
  plotX = 1:length(x);
  %   plotX = (x > -.5) & (x < .5);
  
  tableRowIndex = find(strcmp(psyData.vetmindData{:,1}, target));
  
  for channelIndex = 1:endIndex
    allEvents(:, channelIndex, subjectCounter) = filteredEvent(plotX,channelIndex);
    meanValue = mean(filteredEvent(plotX,channelIndex));
    psyData.vetmindData{tableRowIndex, originalColumnCount + channelIndex} = meanValue;
    if(subjectCounter == 1)
      psyData.vetmindData.Properties.VariableNames{originalColumnCount + channelIndex} = sprintf('%s_phase_height', data.topography.chanlocs(channelIndex).labels);
    end
    meanEvents(channelIndex, subjectCounter) = meanValue;
    if(doPlot)
      close all;
      figure;
      plot(x(plotX), filteredEvent(plotX,channelIndex));
      legend(data.topography.chanlocs(channelIndex).labels);
      title(target)
    end
  end
  
end


%correlations
if(doNeuropsychCorrelations)
  columnCount = size(psyData.vetmindData, 2);
  rhos = NaN(columnCount);
  ps = NaN(columnCount);
  hitCounter = 1;
  for i = 1:columnCount
    for j = (i+1):columnCount
      x = psyData.vetmindData{:, i};
      y = psyData.vetmindData{:, j};
      if(isnumeric(x) && isnumeric(y))
        keep = ~isnan(x) & ~isnan(y);
        [rho, p] = corr(x(keep), y(keep));
        rhos(i,j) = rho;
        rhos(j,i) = rho;
        ps(i,j) = p;
        ps(j,i) = p;
        if(p < 0.05)
          hitP(hitCounter, 1) = p;
          hitRho(hitCounter, 1) = rho;
          hitVar1(hitCounter, 1) = psyData.vetmindData.Properties.VariableNames(i);
          hitVar2(hitCounter, 1) = psyData.vetmindData.Properties.VariableNames(j);
          hitCounter = hitCounter + 1;
        end
      end
    end
  end
  
  tab = table(hitVar1, hitVar2, hitRho, hitP)
  
  figure;
  hit = ps < 0.05;
  imagesc(hit(:, end-31:end));
  
  isFz1 = cellfun(@length, strfind(hitVar1, 'Fz'));
  isFz2 = cellfun(@length, strfind(hitVar2, 'Fz'));
  isFz = isFz1 | isFz2;
  tab(isFz, :)
end


%Interestingly, credibility and expectancy are negatively correlated with
%Fz phase slope.

channelIndex = 0;

channelIndex = channelIndex +1;
toPlot = squeeze(eventData(:, channelIndex, :));
srate = 1024;
x = (((-1*srate):(2.5*srate))./srate)';
% close all;
figure;
plot(x, toPlot);
myTitle= sprintf('%s window size %d', data.topography.chanlocs(channelIndex).labels, data.topography.details.windowSize);
title(myTitle);

%debug
if(data.topography.details.windowSize == 1024)
  largeWindow = eventData;
  save('/home/data/EEG/processed/Oregon/phaseLargeWindow.mat', 'largeWindow', '-v7.3');
else
  smallWindow = eventData;
  save('/home/data/EEG/processed/Oregon/phaseSmallWindow.mat', 'smallWindow', '-v7.3');
end
load('/home/data/EEG/processed/Oregon/phaseLargeWindow.mat');
load('/home/data/EEG/processed/Oregon/phaseSmallWindow.mat');

trialNumber = 0;

trialNumber = trialNumber + 1;
y1 = largeWindow(:,1,trialNumber);
y2 = smallWindow(:,1,trialNumber);
a = 1;
filtSize = 41;
b = repmat(1/filtSize, [1, filtSize]);
y1 = filter(b, a, y1);
y2 = filter(b, a, y2);
y1 = y1 / nanstd(y1);
y2 = y2 / nanstd(y2);

figure;
plot([y1, y2]);
legend({'1 second', '0.05 seconds'});
title(sprintf('trial %d', trialNumber));





