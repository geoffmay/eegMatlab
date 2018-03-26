load('/home/data/EEG/processed/Oregon/reliability3/PM102.mat');

thisDiff = summary.icaCohDiffMeans(:,1);
x = summary.sampleFrameDuration ./ (128 * 60);

figure;
plot(thisDiff);


