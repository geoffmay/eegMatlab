clear;

plotTimeCourse = false;
fullMatrixPlot = false;

folder = '/home/data/EEG/processed/Robi/coherence/';
fileLabels = [{'pre'}, {'post'}];
%folder = '/home/data/EEG/processed/Robi/coherence/ica';
files = dir(folder);
files([files.isdir]) = [];
files=files(find(cellfun(@length, strfind({files.name}, '.mat'))));



for fileCounter=1:length(files)
  filename = fullfile(folder,files(fileCounter).name);
  data = load(filename);
  coherencePlotLength = size(data.channelPairs(1).coherence, 1);
  numberOfChannelPairs = length(data.channelPairs);
  maxChannels = channelCountFromPairs(numberOfChannelPairs);
  
  
  meanSum = zeros(1,5);
  stat.meanCoherences = NaN(maxChannels,maxChannels,5);
  stat.stdDevCoherences = NaN(maxChannels,maxChannels,5);
  stat.interfrequencyCoherenceCorrelation = NaN(maxChannels,maxChannels);
  stat.interfrequencyCoherenceCorrelationP = NaN(maxChannels,maxChannels);
  
  [labels, chanlocs] = antChannelLocs;
  
  channelPairCounter = 1;
  for i = 1:maxChannels
    for j = i+1:maxChannels
      if(length(strfind(files(fileCounter).name, 'ica')) == 0)
        %make sure we're getting what we think we're getting.
        testLabel = sprintf('%s-%s', labels{i}, labels{j});
        test2Label = data.channelPairs(channelPairCounter).label;
        if(~strcmp(testLabel, test2Label))
          error('label mismatch');
        end
      end
      
      coh = data.channelPairs(channelPairCounter).coherence;
      if(plotTimeCourse)
        %debug
        close all;
        figure;
        hold on;
        x = (1:size(coh,1)) ./ 128;
        plot(x, coh);
        zoom xon;
        pan xon;
        title(sprintf('coherence over time, %s-%s', labels{i}, labels{j}));
        legend('delta', 'theta', 'alpha', 'beta', 'high beta');
        xlabel('time (seconds)');
        ylabel('percent phase coherence');
        
        %end debug
      end
      meanCoherences(j,i,:) = mean(coh, 1);
      for k = 1:5
      meanSum(k) = meanSum(k) + meanCoherences(j,i,k);
      end
      if(fullMatrixPlot)
        meanCoherences(i,j,:) = mean(coh, 1);
        meanCoherences(i,i,:) = 1;
        meanCoherences(j,j,:) = 1;
      end
      stdDevCoherences(i,j,:) = std(coh, 1);
      
      rhoSum = 0;
      pSum = 0;
      corrCounter = 0;
      for k = 1:size(coh,2)
        for m = k+1:size(coh,2)
          [rho, p]=corr(coh(:,k),coh(:,m));
          rhoSum = rhoSum + rho;
          pSum = pSum + p;
          corrCounter = corrCounter + 1;
        end
      end
      interfrequencyCoherenceCorrelation(i,j) = rhoSum / corrCounter;
      interfrequencyCoherenceCorrelationP(i,j) = pSum / corrCounter;
      
      channelPairCounter = channelPairCounter + 1;
    end
  end
  stats(fileCounter) = stat;
  if(false)
    close all;
    for mapNumber = 1:size(data.ica.mixingMatrix, 2)
      figure;
      map = data.ica.mixingMatrix(:,mapNumber);
      topoplot(map,chanlocs(1:32))
      title(sprintf('component %d', mapNumber));
    end
    tilefigs;
  end
  %color unimportant boxes with the mean
  if(~fullMatrixPlot)
    for i = 1:maxChannels
      for j = i:maxChannels
        for k = 1:5
        meanCoherences(i,j,k) = meanSum(k)./numberOfChannelPairs;
        end
      end
    end
  end
  
  close all;
  figure;
  imagesc(meanCoherences(:,:,5));
  title(sprintf('beta coherence %s', fileLabels{fileCounter}));
  set(gca, 'xtick', 1:32);
  set(gca, 'xticklabel', labels(1:32));
  set(gca, 'ytick', 1:32);
  set(gca, 'yticklabel', labels(1:32));
  handleToColorbar = colorbar;
  title(handleToColorbar, 'phase coherence');
  
end




% for i = 1:length(data.channelPairs)
%   coh = data.channelPairs(i).coherence;
%   betaCoh(:,i) = coh(:,5);
% %   figure;
% %   plot(betaCoh(:,i));
% %   title(pairLabels{i});
% end
% remove = all(isnan(betaCoh));
% betaCoh(:,remove) = [];
% labels(remove) = [];
%
%
%
% remove = all(isnan(betaCoh),2);
% betaCoh(find(remove),:) = [];
%




%
% %outcome with different structs
% if(false)
% betaCoh = NaN(size(coherences{1},1),length(coherences));
%
% for i = 1:length(coherences)
%   coh = coherences{i};
%   betaCoh(:,i) = coh(:,5);
% %   figure;
% %   plot(betaCoh(:,i));
% %   title(pairLabels{i});
% end
% labels = pairLabels;
%
% remove = all(isnan(betaCoh));
% betaCoh(:,remove) = [];
% labels(remove) = [];
%
%
%
% remove = all(isnan(betaCoh),2);
% betaCoh(find(remove),:) = [];
% end
%
% %dendrogram stuff
% if(false)
% tic;
% z1 = linkage(betaCoh', 'average', 'euclidean');
% toc;
% [H, T] = dendrogram(z1, 0);
% c = garlicCluster(z1, .005, .2);
% sortC = sortrows(tabulate(c), 2);
% biggestGroup = sortC(end,1);
% indices = find(c == biggestGroup);
% indices1 = T(indices);
% clusterLabels = cell(0);
% for i = 1:length(indices1)
%   index1 = find(z1(:,1)==indices1(i));
%   index2 = find(z1(:,2)==indices1(i));
%   if(length(index1) > 0)
%     if(index1 <= length(labels))
%       clusterLabels{end+1} = labels{index1};
%     end
%   end
%   if(length(index2) > 0)
%     if(length(index2) > 0)
%       clusterLabels{end+1} = labels{index2};
%     end
%   end
% end
% uniqueLabels = unique(clusterLabels);
%
% end
%
% %huge covariance matrix
% if(false)
% close all;
% pairBetaCovariance = cov(betaCoh);
% imagesc(pairBetaCovariance);
% end