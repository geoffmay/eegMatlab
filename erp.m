function [output] = erp()

allEvents = [];
load(fullfile('/home/gmay/Documents/MATLAB', 'eventData.mat'));

tabulation = tabulate([allEvents.commonEvent]);
tabulation(find(tabulation(:,2) == 0),:) = [];
tabulation = sortrows(tabulation,2);
mostCommonEvent = tabulation(end,1);
mostCommonEventCount = tabulation(end,2);

eventsOfInterest = [129:132, 4];

todo = find(cellfun(@length, strfind({allEvents.filename}, 'Flanker')));

%for i = 1:length(allEvents)
for fileIndex = 1:length(todo)
  i = todo(fileIndex);
  %   display(strcat('processing file (', num2str(i),...
  %     ') of (', num2str(mostCommonEventCount), ')'));
  a = allEvents(i);
  %    if(a.commonEvent == mostCommonEvent)
  slashes = strfind(allEvents(i).filename, '/');
  file = allEvents(i).filename(slashes(end-1):end);
  folder = '/home/data/EEG/data/Oregon';
  path = fullfile(folder, file);
  data.file = path;
     [sourceFolder, sourceFile] = fileparts(file);
    %if(~exist(destination))       
    try
      EEG = pop_readbdf(path, {}, 43, int32(32), false);
      for eventCounter = 1:length(eventsOfInterest);
        data.event = eventsOfInterest(eventCounter);
        hits = find(a.sequence(:,1) == data.event);
        if(length(hits) > 0)
          
          for j = 1:length(EEG.chanlocs)
            fprintf('file %d of %d, channel %d of %d\n', fileIndex, length(todo), j, length(EEG.chanlocs));
            %         channel.fullTimeF = waveletAnalysis(EEG, [2, 0.5], data.event, EEG.chanlocs(j).labels);
            
            channel.label = EEG.chanlocs(j).labels;
            data.channels(j) = channel;
          end          
data
summaries(eventCounter) = data;
        end        
      end
      destFile = sprintf('%s%s%s', sourceFile, 'event', '.mat');
      destination = sprintf('/home/data/EEG/processed/Oregon/ERP/summary%d', data.event);
      if(~exist(destination))
        mkdir(destination);
      end
      destination = fullfile(destination, destFile);
      save(destination, 'summaries');
    catch
      disp('error processing file');
    end
    %end
  %data.timeF = waveletAnalysis(EEG, [2, 0.5], 209, 'Fz');
  %         save(fullfile(folder, 'E:\EEG\timeFreqDecomp.mat', 'fullRecord'));
  %         picName = strcat('E:\EEG\timeF\', files(i).filename, 'TimeF.png');
  %         print(picName, '-dpng');
  fclose('all');
  %         close all;
  %    end
end

folder = '/home/data/EEG/data/wavelet';
outputFolder = '/home/data/EEG/processed/fullTimeF';
fileNames = dir(folder);
for i = length(fileNames):-1:1
  if(fileNames(i).name(1) == '.' || ~strcmp(fileNames(i).name(end-3:end), '.mat'))
    fileNames(i) = [];
  end
end

channelLabels = wahbehChannelLocs;

searchString = 'AF3.mat';
batchFilenameIndent = strfind({fileNames.name}, searchString);
batchFileNumber = find(cellfun(@length, batchFilenameIndent) > 0);
batchFileNumber(end+1) = length(fileNames) + 1;
batchSpacing = diff(batchFileNumber);
batchFileNumber(end) = [];
channelCount = mode(batchSpacing);
%files = cell(1,channelCount);
for i = 1:length(batchFileNumber)
  i1 = batchFileNumber(i);
  if(batchSpacing(i) == channelCount) %only process complete sets
    prefix = fileNames(i1).name(1:batchFilenameIndent{i1}-1);
    outputFilename = fullfile(outputFolder,strcat(prefix,'.csv'));
    if(~exist(outputFilename))
      disp(sprintf('%d of %d', i, length(batchFileNumber)));
      for j = 1:length(channelLabels)
        filename = fullfile(folder, strcat(prefix, channelLabels{j}, '.mat'));
        if(~exist(filename))
          error('no file.');
        else
          load(filename, 'fullTimeF');
          if(j == 1)
            allErspDims = size(fullTimeF.allErsp{1});
            width = allErspDims(1) * allErspDims(2) * channelCount;
            data = NaN(length(fullTimeF.allErsp), width);
            
          end
          for k = 1:length(fullTimeF.allErsp)
            m = fullTimeF.allErsp{k};
            mf = reshape(m, 1, size(m,1)*size(m,2));
            x = (j-1) * length(mf) + 1;
            data(k, x:x+length(mf)-1) = mf;
          end
        end
      end
      csvwrite(outputFilename, data);
    end
  end
end




end