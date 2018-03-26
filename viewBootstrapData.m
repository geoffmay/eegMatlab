folder = '/home/data/EEG/processed/Robi/coherenceBootstrap';
filenames = dir(folder);
filenames([filenames.isdir]) = [];
for fileCounter = 1:length(filenames)
  filename = filenames(fileCounter).name;
  load(filename);
  for i = 1:length(data.ZscoreFile1)
    maxZ1(i) = max(abs(data.ZscoreFile1(i)));
    maxZ2(i) = max(abs(data.ZscoreFile2(i)));
    p(i) = normcdf(min(-maxZ1(i), -maxZ2(i)));
    logp(i) = log10(p(i));
  end
  fprintf('\n%f: %s', max(logp), filename);
  
  
  for i = 1:length(data.channelSummaries)
    a1 = data.channelSummaries(i).Zscore1(end);
    a2 = data.channelSummaries(i).Zscore2(end);
    z(i) = max(abs(a1),abs(a2));
    chanP(i) = normcdf(-z(i));
    logChanP(i) = log(chanP(i));
  end
  
  if(false)
    baseline = load('ROBI_003_baseline eyes open_630158995692243270coherenceStats');
    outcome = load('ROBI_003_outcome eyes open_630230539149591228coherenceStats');
    [labels,chanlocs] = antChannelLocs;
    m1 = find(strcmp(labels, 'M1'));
    m2 = find(strcmp(labels, 'M2'));
    remove = [m1 m2];
    baseline.stat.meanCoherences(remove,:,:) = [];
    baseline.stat.meanCoherences(:,remove,:) = [];
    outcome.stat.meanCoherences(remove,:,:) = [];
    outcome.stat.meanCoherences(:,remove,:) = [];
    labels(remove) = [];
    chanlocs(remove) = [];
    difference = outcome.stat.meanCoherences - baseline.stat.meanCoherences;
    counter = zeros(1, size(difference,1));
    sums = zeros(1, size(difference,1));
    for i = 1:size(difference,1)
      for j = (i+1):size(difference,2)
        %for k = 1:size(difference,3)
        for k = 1:5
          difference(j,i,k) = difference(i,j,k);
          sums(i) = sums(i) + difference(i,j,k);
          sums(j) = sums(j) + difference(i,j,k);
          counter(i) = counter(i) + 1;
          counter(j) = counter(j) + 1;
        end
      end
    end
    
    for i = 1:length(sums)
      avgs(i) = sums(i) / counter(i);
    end
    figure;
    topoplot(avgs.*100, chanlocs(1:length(avgs)), 'electrodes', 'labels','gridscale',300);
    title('global change in coherence, pre- to post-neurofeedback');
    cbar = colorbar;
    title(cbar, 'absolute change (%)');
    meanChange = mean(mean(difference, 1), 3);
    
    meanChange2 = squeeze(meanChange);
  end
  
  for i = 1:32
  end
  
  [labels, chanlocs] = antChannelLocs;
  
  freqLabels = [{'delta'}, {'theta'}, {'alpha'}, {'beta'}, {'hibeta'}, {'gamma'}, {'all'}];
  tab = table(freqLabels', logp');
  
  correctionFactor = 6 * 31 * (31-1) / 2;
  bonferroniThreshold = log10(0.05 / correctionFactor);
  
  percentiles = NaN(33,1);
  freqLabels = [{'delta'},{'theta'},{'alpha'},{'beta'},{'hibeta'},{'gamma'},{'average'}];
  frequencyBand = 4;
  for i = 1:length(data.channelSummaries)
    p = data.channelSummaries(i).file1Percentile(frequencyBand);
    inverter = 1;
    if(p > .5)
      inverter = -1;
      p = 1 - p;
    end
    p = log10(p);
    if(isinf(p))
      p = -6;
    end
    if(p > bonferroniThreshold)
      p = 0;
    end
    p = p * inverter;
    percentiles(i) = p;
  end
  m1 = find(strcmp(labels,'M1'));
  m2 = find(strcmp(labels,'M2'));
  remove = [m1 m2];
  percentiles(remove) = [];
  chanlocs(remove) = [];
  chanlocs = chanlocs(1:length(percentiles));
  figure;
  topoplot(percentiles, chanlocs, 'electrodes', 'labels');
  subjectName = filename(1:18);
  subjectName = strrep(subjectName,'_',' ');
  title(sprintf('%s coherence changes %s', freqLabels{frequencyBand}, subjectName));
  caxis([-6 6]);
  colorbar;
  
  chanPairs = zeros(33,33);
  i = 0;
  for chan1 = 1:33
    for chan2 = (chan1+1):33
      i = i + 1;
      p = data.channelPairs(i).file1Percentile(frequencyBand);
      inverter = 1;
      if(p > .5)
        inverter = -1;
        p = 1 - p;
      end
      p = log10(p);
      if(isinf(p))
        p = -6;
      end
      p = p * inverter;
      chanPairs(chan1,chan2) = p;
      chanPairs(chan2,chan1) = p;
    end
  end
  
  chanPairs1 = chanPairs;
  chanPairs1(remove,:) = [];
  chanPairs1(:,remove) = [];
  
  
  figure;
  imagesc(chanPairs1);
  cmin = min(percentiles);
  cmax = max(percentiles);
  if(cmin==cmax)
    cmax = cmax + 1;
  end
%  caxis([cmin, cmax]);
  caxis([-6 6]);
  cbar = colorbar();
  title(cbar, 'log10p');
  title(sprintf('%s coherence changes from pre to post', freqLabels{frequencyBand}));
  set(gca, 'xtick', 1:33);
  set(gca, 'xticklabel', {chanlocs.labels});
  set(gca, 'ytick', 1:33);
  set(gca, 'yticklabel', {chanlocs.labels});
  
end
