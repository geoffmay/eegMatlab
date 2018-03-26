clear;

overwrite = false;
% filenames = [...
%     {'E:\EEG\ROBI\raw data\ROBI_003\baseline eyes open\630158995692243270.eegData'},...
%     {'E:\EEG\ROBI\raw data\ROBI_003\outcome eyes open\630230539149591228.eegData'}...
%     ];

filenames = getRobiDataFiles();
timeStamps = NaN();
for i = 1:length(filenames)
  filename = filenames{i};
  slashes = strfind(filenames{i}, '/');
  numberEnd = strfind(filenames{i}, '.eegData')-1;
  numberStart = slashes(end)+1;
  subString = filename(numberStart:numberEnd);
  timeStamps(i) = str2num(subString);
end
fileTable = table(filenames', timeStamps');
sortedFiles = sortrows(fileTable,2);

channelCount = 34;


for fileCounter = 1:length(filenames)
  filename = sortedFiles{fileCounter,1};
  filename = filename{1};
  outputFolder = '/home/data/EEG/processed/Robi/phaseSlope/';
  outputStart = strfind(filename, 'ROBI_');
  if(length(outputStart) < 1)
    error('invalid filename');
  end
  outFile = filename(outputStart(1):end);
  outFile = strrep(outFile,'/','_');
  outFile = strrep(outFile,'.eegData','phaseSlope.mat');
  outputPath = fullfile(outputFolder, outFile);
  
  write = false;
  if(~exist(outputPath,'file'))
    write = true;
  elseif(overwrite)
    write = true;
  end
  
  if(write)    
    %Info.mat';
    fprintf('%s', filename);
    file = dir(filename);
    fileLength = file.bytes / 8;
    summary.filename = filename;
    
    fileId = fopen(filename, 'r');
    contents = fread(fileId, fileLength, 'double');
    sampleCount = fileLength / channelCount;    
    if(sampleCount ~= floor(sampleCount))
      sampleCount = floor(sampleCount);
      fileLength = sampleCount * channelCount;
      contents = contents(1:fileLength);
    end

    fclose(fileId);    
    data = reshape(contents, channelCount, fileLength / channelCount)';
    channelCounter = 1;
    labels = antChannelLocs;
    
    for chanIndex1 = 1:32
      for chanIndex2 = chanIndex1+1:32
        %             fprintf('.');
        %             if(mod(channelCounter,100) == 0)
        %                 fprintf('\n%d',size(channelPair,1));
        %             end
        channelPair.channel1 = labels{chanIndex1};
        channelPair.channel2 = labels{chanIndex2};
        tic;
        [channelPair.slope, channelPair.phaseAngle, channelPair.poly] = phaseSlopeIndex(data, chanIndex1, chanIndex2);
        elapsed = toc;
        fprintf('\n%d/%d(%d)%s-%s: %f [%f seconds]', fileCounter,length(filenames), channelCounter, channelPair.channel1, channelPair.channel2, channelPair.slope, elapsed);
        if(~exist('channelPairs', 'var'))
          channelPairs = channelPair;
        else
          channelPairs(end+1) = channelPair;
        end
        %debug
        close all;
        figure;
        hold on;
        plot(channelPair.phaseAngle);
        polyPredict = (2:40) .* channelPair.poly(1) + channelPair.poly(2);
        plot(2:40, polyPredict, 'r');
        title(sprintf('\n(%d)%s-%s: %f', channelCounter, channelPair.channel1, channelPair.channel2, channelPair.slope));
        %end debug
        channelCounter = channelCounter + 1;
      end
    end
    
    
    summary.channelPairs = channelPairs;
    save(outputPath, 'summary', '-v7.3');
    clear channelPairs;
    %     if(~exist('summaries', 'var'))
    %         summaries = summary;
    %     else
    %         summaries(end+1) = summary;
    %     end
    clear summary;
  end
end
