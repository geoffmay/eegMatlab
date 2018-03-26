cd /home/data/EEG/scripts/bespoke/;
outputFolder = '/home/data/EEG/processed/Oregon/reliability2';

snipDurations = 20:20:120;
repetitions = 3;
fileName = '/home/data/EEG/data/Oregon/PTSD/PM101.bdf';
[folder, file, ext] = fileparts(fileName);
outputPath = fullfile(outputFolder, sprintf('%sSnippets.mat', file));
eeg = loadBdf(fileName);
snip = eeg;

%loop through different snip durations
clear summaries;
rhos = NaN(repetitions, length(snipDurations));
for snipDurIndex = 1:length(snipDurations)
  sampleDuration = snipDurations(snipDurIndex);
  maxStart = size(eeg.data, 2) - sampleDuration * eeg.srate;
  start = ceil(rand * maxStart);
  snip.data = eeg.data(:, start:(start+sampleDuration*eeg.srate-1));
  for i = 1:repetitions
    timer = tic;
    summary = performReliabilityAnalysis(snip);
    summary.snippetDuration = sampleDuration;
    summary.sampleRate = 128;
    summary.secondsNeededToCompute = toc(timer);
    if(false)
      meanRho = mean(summary.icaRhos, 2);
      meanSurfRho = mean(summary.surfRhos, 2);
      x = summary.sampleFrameDuration ./ 128;
      figure;
      plot(x, [meanRho meanSurfRho])
      legend('ica', 'surface');
      xlabel('sub-sample duration (seconds)');
      ylabel('Pearson''s rho');
      title(sprintf('reliability of sub-samples in a %d second eeg clip', sampleDuration));
      
      figure;
      plot(x, summary.icaRhos);
    end
    summaries(i, snipDurIndex) = summary;
  end  
end