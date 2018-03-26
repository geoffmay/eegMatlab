function data = bootstrap(file1, file2)

bootstrapSampleSize = 10000;
bootstrapRepetitionCount = 100000;
percentiles = [0.0001, 0.001, 0.05, 0.95, 0.999, 0.9999];

disp('loading...');
if(~exist('file1','var'))
  file1 = '/home/data/EEG/data/ROBI/ROBI_003/baseline eyes open/630158995692243270.eegData';
  file2 = '/home/data/EEG/data/ROBI/ROBI_003/outcome eyes open/630230539149591228.eegData';
end
data1 = loadRobiDataFile(file1);
data2 = loadRobiDataFile(file2);
data.file1 = file1;
data.file2 = file2;

data.lowFreq = [1 4 8 12 25 30];
data.highFreq = [4 8 12 25 30 60];


labels = antChannelLocs;
pairCounter = 1;
% channelSums1 = zeros(2,2,33);
% channelSums2 = zeros(2,2,33);
grandSum1 = 0;
grandSum2 = 0;
for chan1 = 1:33
  for chan2 = (chan1+1):33
    fprintf('\n%s, pair %d',char(datetime), pairCounter);
    %compute coherence for the channel pair
    coh1 = coherence(data1(:,chan1), data1(:,chan2), data.lowFreq, data.highFreq);
    coh2 = coherence(data2(:,chan1), data2(:,chan2), data.lowFreq, data.highFreq);
    coh1 = [coh1 mean(coh1,2)];
    coh2 = [coh2 mean(coh2,2)];
    
    %truncate to permit maximum possible variance for both samples
    if(size(coh1,1) > size(coh2,1))
      coh1((size(coh2,1)+1):end,:) = [];
    elseif(size(coh1,1) < size(coh2,1))
      coh2((size(coh1,1)+1):end,:) = [];
    end
    
    %downsample to remove autocorrelation
    downSize = floor(size(coh1,1) / 512);
    tempC1 = NaN(downSize, size(coh1,2));
    tempC2 = NaN(downSize, size(coh2,2));
    sampler = 1;
    for downSampleCounter = 1:downSize
      tempC1(downSampleCounter,:) = coh1(sampler, :);
      tempC2(downSampleCounter,:) = coh2(sampler, :);
      sampler = sampler + 512;
    end
    coh1 = tempC1;
    coh2 = tempC2;
    
    %update sums
    if(chan1 == 1)
      if(chan2 == 2)
        channelSums1 = zeros(size(coh1,1),size(coh1,2),33);
        channelSums2 = zeros(size(coh1,1),size(coh1,2),33);
        channelSums1(:,:,chan1) = coh1;
        channelSums2(:,:,chan1) = coh2;
      else
        channelSums1(:,:,chan1) = channelSums1(:,:,chan1) + coh1;
        channelSums2(:,:,chan1) = channelSums2(:,:,chan1) + coh2;
      end
      channelSums1(:,:,chan2) = coh1;
      channelSums2(:,:,chan2) = coh2;
    else
      channelSums1(:,:,chan1) = channelSums1(:,:,chan1) + coh1;
      channelSums2(:,:,chan1) = channelSums2(:,:,chan1) + coh2;
      channelSums1(:,:,chan2) = channelSums1(:,:,chan2) + coh1;
      channelSums2(:,:,chan2) = channelSums2(:,:,chan2) + coh2;
    end
    grandSum1 = grandSum1 + coh1;
    grandSum2 = grandSum2 + coh2;

    %bootstrap individual channel pair
    combined = [coh1; coh2];
    sampleStats = NaN(bootstrapRepetitionCount, size(combined,2));
    for repetitionCounter = 1:bootstrapRepetitionCount
      sampleNumbers = ceil(rand(1,bootstrapSampleSize) .* size(combined,1));
      resampledData = combined(sampleNumbers,:);
      sampleStats(repetitionCounter,:) = mean(resampledData,1);
    end
    
    %summarize percentiles
    a = array2table([percentiles', prctile(sampleStats, percentiles, 1)]);
    a.Properties.VariableNames{1} = 'percentile';
    a.Properties.VariableNames{2} = 'delta';
    a.Properties.VariableNames{3} = 'theta';
    a.Properties.VariableNames{4} = 'alpha';
    a.Properties.VariableNames{5} = 'beta';
    a.Properties.VariableNames{6} = 'hibeta';
    a.Properties.VariableNames{7} = 'gamma';
    a.Properties.VariableNames{8} = 'allFrequencies';
    channelPair.percentiles = a;
    
    %compute percentile for file1 and file2.
    cohMean1 = mean(coh1,1);
    cohMean2 = mean(coh2,1);
    for i = 1:size(cohMean1,2)
      less1 = cohMean1(i) < sampleStats(:,i);
      channelPair.file1Percentile(i) = sum(less1)/length(less1);
      less2 = cohMean2(i) < sampleStats(:,i);
      channelPair.file2Percentile(i) = sum(less2)/length(less2);      
      %addition begins here
      sampleMean = mean(sampleStats(:,i));
      sampleStd = std(sampleStats(:,i));
      channelPair.bootstrapSampleMean(i) = sampleMean;
      channelPair.bootstrapSampleStd(i) = sampleStd;
      channelPair.Zscore1(i) = (cohMean1(i) - sampleMean) / sampleStd;
      channelPair.Zscore2(i) = (cohMean2(i) - sampleMean) / sampleStd;
    end
    
    channelPair.chan1 = labels{chan1};
    channelPair.chan2 = labels{chan2};
    data.channelPairs(pairCounter) = channelPair;
    pairCounter = pairCounter + 1;
  end
  
  fprintf('\n%s, channel %d',char(datetime), chan1);
  %bootstrap channel summary
  channelSummary.label = labels{chan1};
  sum1 = channelSums1(:,:,chan1);
  sum2 = channelSums2(:,:,chan1);
  combined = [sum1;sum2];
  sampleStats = NaN(bootstrapRepetitionCount, size(combined,2));
  for repetitionCounter = 1:bootstrapRepetitionCount
    sampleNumbers = ceil(rand(1,bootstrapSampleSize) .* size(combined,1));
    resampledData = combined(sampleNumbers,:);
    sampleStats(repetitionCounter,:) = mean(resampledData,1);
  end
  cohMean1 = mean(sum1,1);
  cohMean2 = mean(sum2,1);
  for i = 1:size(cohMean1,2)
    less1 = cohMean1(i) < sampleStats(:,i);
    channelSummary.file1Percentile(i) = sum(less1)/length(less1);
    less2 = cohMean2(i) < sampleStats(:,i);
    channelSummary.file2Percentile(i) = sum(less2)/length(less2);
    %addition begins here
    sampleMean = mean(sampleStats(:,i));
    sampleStd = std(sampleStats(:,i));
    channelSummary.bootstrapSampleMean(i) = sampleMean;
    channelSummary.bootstrapSampleStd(i) = sampleStd;
    channelSummary.Zscore1(i) = (cohMean1(i) - sampleMean) / sampleStd;
    channelSummary.Zscore2(i) = (cohMean2(i) - sampleMean) / sampleStd;
  end
  data.channelSummaries(chan1) = channelSummary;  
end
fprintf('\n%s, mean bootstrap',char(datetime));
%bootstrap grand mean
combined = [grandSum1;grandSum2];
sampleStats = NaN(bootstrapRepetitionCount, size(combined,2));
for repetitionCounter = 1:bootstrapRepetitionCount
  sampleNumbers = ceil(rand(1,bootstrapSampleSize) .* size(combined,1));
  resampledData = combined(sampleNumbers,:);
  sampleStats(repetitionCounter,:) = mean(resampledData,1);
end
cohMean1 = mean(grandSum1,1);
cohMean2 = mean(grandSum2,1);
for i = 1:size(cohMean1,2)
  less1 = cohMean1(i) < sampleStats(:,i);
  data.file1Percentile(i) = sum(less1)/length(less1);
  less2 = cohMean2(i) < sampleStats(:,i);
  data.file2Percentile(i) = sum(less2)/length(less2);
  %addition begins here
  sampleMean = mean(sampleStats(:,i));
  sampleStd = std(sampleStats(:,i));
  data.bootstrapSampleMean(i) = sampleMean;
  data.bootstrapSampleStd(i) = sampleStd;
  data.ZscoreFile1(i) = (cohMean1(i) - sampleMean) / sampleStd;
  data.ZscoreFile2(i) = (cohMean2(i) - sampleMean) / sampleStd;
end

end