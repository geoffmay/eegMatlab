folder = '/home/data/EEG/processed/Robi/coherence/stats';
files = dir(folder);
files([files.isdir]) = [];
counter = 1;
for i = 1:length(files)
  
  filename = files(i).name;
  data = load(fullfile(folder,filename));
  if(isfield(data, 'stat'))
    datas(counter,1) = data.stat;
    underscores = strfind(filename,'_');
    startNumber = underscores(end) + 1;
    endNumber = startNumber + length('630182103318938594')-1;
    numberString = filename(startNumber:endNumber);
    numbers(counter,1) = str2num(numberString);
    filenames{counter,1} = data.stat.filename;
    counter = counter + 1;
  end
end
tab = table(datas, numbers, filenames);
filter = 'ROBI_003';
filterMiss = find(~cellfun(@length,strfind(tab{:,'filenames'}, filter)));
tab(filterMiss,:) = [];
tab = sortrows(tab,2);
tab(:,2:3)=[];
for fileCounter = 1:size(tab,1)
  data = tab{fileCounter,1};
  meanCoh = data.meanCoherences;
  meanCoh(isnan(meanCoh)) = 0;
  meanmean = mean(mean(meanCoh));
  grandMean = mean(meanmean);
  deltaMean(fileCounter,1) = meanmean(:,:,1);
  thetaMean(fileCounter,1) = meanmean(:,:,2);
  alphaMean(fileCounter,1) = meanmean(:,:,3);
  betaMean(fileCounter,1) = meanmean(:,:,4);
  hibetaMean(fileCounter,1) = meanmean(:,:,5);
  
  meanCoherence(fileCounter,1) = grandMean;
end

tab2 = table(meanCoherence, deltaMean, thetaMean, alphaMean, betaMean, hibetaMean);

close all;
figure;
hold on;
plot(tab2{:,:});
labels = tab2.Properties.VariableNames;
legend(labels);