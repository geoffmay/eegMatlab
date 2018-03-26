clear;

createDummyFile = true;

doPhaseAngle = false;
simOnly = false;
summaryOnly = false;
removeMean = true;


sampleRate = 2048;
epochLength = 1 * sampleRate;
chan1 = 100;
startIndex = (chan1-1) * epochLength + 1;
endIndex = startIndex + epochLength - 1;

hiPassHz = 1;
loPassHz = 40;

% filenames = [{'/home/data/EEG/data/ROBI/ROBI_003/baseline eyes open/630158995692243270.eegData'}, ...
%   {'/home/data/EEG/data/ROBI/ROBI_003/outcome eyes open/630230539149591228.eegData'}];

%filename = '/Users/Geoff/Documents/MATLAB/EEG/Coe Collection/Robi/ROBI_003/tx 1/630165006453007385.eegData';
%        filename = '/home/data/EEG/data/ROBI/ROBI_003/outcome eyes open/630230539149591228.eegData';

filenames = getRobiDataFiles();
lowFreq = [1 4 8 12 25];
highFreq = [4 8 12 25 30];

for fileCounter = 1:length(filenames)
  filename = filenames{fileCounter};
  if(summaryOnly)
  outputFolder = '/home/data/EEG/processed/Robi/power';
  else
      outputFolder = '/media/eegDrive/power';
  end
  outputStart = strfind(filename, 'ROBI_');
  if(length(outputStart) < 1)
    error('invalid filename');
  end
  outFile = filename(outputStart(1):end);
  outFile = strrep(outFile,'/','_');
  outFile = strrep(outFile,'.eegData','PowerStats.mat');
  outputPath = fullfile(outputFolder, outFile);
  %if(~exist(outputPath, 'file'))
    if(createDummyFile)
      dummy = sprintf('started on %s', char(datetime));
      save(outputPath, 'dummy');
    end
    
    channelCount = 34;
    file = dir(filename);
    fileLength = file.bytes / 8;
    fileId = fopen(filename);    
    contents = fread(fileId, fileLength, 'double');    
    fclose(fileId);
    
    %truncate if necessary
    sampleCount = fileLength / channelCount;
    if(sampleCount ~= floor(sampleCount))
      sampleCount = floor(sampleCount);
      fileLength = sampleCount * channelCount;
      contents = contents(1:fileLength);
    end
    
    
    
    labels = antChannelLocs;
    
    data = reshape(contents, channelCount, fileLength / channelCount)';
    clear contents;
    if(removeMean)
        altData = data(:, 1:33);
        meanData = mean(altData, 2);
        for i = 1:size(altData,2)
            altData(:,i) = altData(:,i) - meanData;
        end
        data = altData;
    end
    cpzIndex = find(strcmp(labels,'CPz'));
    m1Index = find(strcmp(labels,'M1'));
    m2Index = find(strcmp(labels,'M2'));
    
    for i=1:size(data,1)
      sample = data(i,:);
      linkedMastoids = (sample(m1Index) + sample(m2Index)) / 2;
      avg = mean(sample(1:33));
      newSample = sample - linkedMastoids;
      avgSample = sample - avg;
      data(i,1:33) = newSample(1:33);      
    end

    
    maxPairs = 0;
    indexTable = [];
    pairLabels = cell(0);
    maxChannel = channelCount - 2;
    for chan1 = 1:(maxChannel)
      for chan2 = chan1+1:(maxChannel)
        maxPairs = maxPairs + 1;
        indexTable(maxPairs, 1) = maxPairs;
        indexTable(maxPairs, 2) = chan1;
        indexTable(maxPairs, 3) = chan2;
        pairLabels{maxPairs}= sprintf('%s-%s',labels{chan1},labels{chan2});
      end
    end
    
    channelPair.asymmetry = [];
    channelPair.label = '';
    pairCounter = 1;
    
    stat.meanAsymmetries = NaN(maxChannel,maxChannel,5);
    stat.stdDevAsymmetries = NaN(maxChannel,maxChannel,5);
    stat.skewAsymmetries = NaN(maxChannel,maxChannel,5);
    stat.kurtosisAsymmetries = NaN(maxChannel,maxChannel,5);

    stat.meanAbsPowers = NaN(maxChannel,5);
    stat.stdDevAbsPowers = NaN(maxChannel,5);
    stat.skewAbsPowers = NaN(maxChannel,5);
    stat.kurtosisAbsPowers = NaN(maxChannel,5);

    stat.meanRelPowers = NaN(maxChannel,5);
    stat.stdDevRelPowers = NaN(maxChannel,5);
    stat.skewRelPowers = NaN(maxChannel,5);
    stat.kurtosisRelPowers = NaN(maxChannel,5);

    stat.meanPowerRatios = NaN(maxChannel,5);
    stat.stdDevPowerRatios = NaN(maxChannel,5);
    stat.skewPowerRatios = NaN(maxChannel,5);
    stat.kurtosisPowerRatios = NaN(maxChannel,5);

    stat.channelLabels = labels(1:32);
    stat.filename = filename;
    
    windowSize = sampleRate;
    
    fprintf('\ncomputing power');
    for chan1 = 1:(maxChannel)
      fprintf('.');
      windowCounter = 1;
      windowStart = 1;
      windowEnd = windowStart + windowSize - 1;
      while(windowEnd < size(data,1))
        window = data(windowStart:windowEnd, chan1);        
        myFft = fft(window);
        sigPower = abs(myFft(1:highFreq(end)+1));
        for(freqCounter = 1:length(lowFreq))
          lo = lowFreq(freqCounter) + 1;
          hi = highFreq(freqCounter) + 1;
          powers(windowCounter, freqCounter) = mean(sigPower(lo:hi));
        end       
        windowStart = windowStart + 16;
        windowEnd = windowStart + windowSize - 1;
        windowCounter = windowCounter + 1;
      end
      if(chan1 == 1)
        chanPowers = NaN(size(powers,1),size(powers,2),maxChannel);
      end
      chanPowers(:,:,chan1) = powers;
      stat.meanAbsPowers(chan1,:) = log10(mean(powers, 1)) + 3;
      stat.stdDevAbsPowers(chan1,:) = log10(std(powers, 1)) + 3;
      stat.skewAbsPowers(chan1,:) = log10(skewness(powers, 1)) + 3;
      stat.kurtosisAbsPowers(chan1,:) = log10(kurtosis(powers, 1)) + 3;
      relPowers = powers;
      for i = 1:size(relPowers,1)
        totalPower = relPowers(i,:);
        relPowers(i,:) = relPowers(i,:) ./ totalPower;
      end
      stat.meanRelPowers(chan1,:) = log10(mean(relPowers,1)) + 3;
      stat.stdDevAbsPowers(chan1,:) = log10(std(relPowers,1)) + 3;
      stat.skewAbsPowers(chan1,:) = log10(skewness(relPowers,1)) + 3;
      stat.kurtosisAbsPowers(chan1,:) = log10(kurtosis(relPowers,1)) + 3;

    end
    
    if(summaryOnly)
    for chan1 = 1:maxChannel
      for chan2 = 1:maxChannel
        numerator = chanPowers(:,:,chan1) - chanPowers(:,:,chan2);
        denominator = chanPowers(:,:,chan1) + chanPowers(:,:,chan2);
        asymmetry = 200 .* numerator ./ denominator;
        stat.meanAsymmetries(chan1, chan2, :) = mean(asymmetry,1);
        stat.stdDevAsymmetries(chan1, chan2, :) = std(asymmetry,1);
        stat.skewAsymmetries(chan1, chan2, :) = skewness(asymmetry,1);
        stat.kurtosisAsymmetries(chan1, chan2, :) = kurtosis(asymmetry,1);
      end
    end
    end
    
    if(summaryOnly)
        save(outputPath, 'stat', '-v7.3');
    else
        save(outputPath, 'chanPowers', '-v7.3');
    end
  %end
end



%
%
%     for pairCounter = 1:maxPairs
%       i = 1;
%       j = 2;
%       thisPair = 1;
%       while(thisPair ~= pairCounter)
%         j = j + 1;
%         if(j > length(labels))
%           i = i + 1;
%           j = i + 1;
%         end
%         thisPair = thisPair + 1;
%       end
%       fprintf('%d, %d\n',i,j);
%
% %       i = indexTable(pairCounter,2);
% %       j = indexTable(pairCounter,3);
%       channelPair(pairCounter).coherence = coherence(data(:,i), data(:,j));
%       channelPair(pairCounter).label = labels(pairCounter);
%       fprintf('%d/%d\n',pairCounter,maxPairs);
%       channelPairs(pairCounter) = channelPair;
%     end
%     save('/home/data/EEG/processed/Robi/coherence/ROBI_003_outcomeEOCoherence.mat', 'channelPairs');
%
