filename = '/home/data/EEG/processed/Robi/coherenceBootstrap/ROBI_003_eyes open bootstrap.mat';
load(filename);

displayValues = NaN(length(data.channelSummaries),2);
for i = 1:size(displayValues,1)
  p(1) = data.channelSummaries(i).file1Percentile;
  p(2) = data.channelSummaries(i).file2Percentile;
  for j = 1:length(p)    
    q = p(j);
    inverter = 1;
    if(q > .5)
      inverter = -1;
      q = 1-q;
    end
    d = log10(q);
    if(isinf(d))
      d = -6;
    end
    d = d * inverter;
    displayValues(i,j) = d;
  end
end

labels = antChannelLocs;
if(~exist('topoplot','file'))
  eeglab;
end
topoplot(displayValues, labels);

