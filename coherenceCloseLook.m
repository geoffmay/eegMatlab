if(~exist('channelPairs', 'var'))
folder = '/media/eegDrive';
filename = 'ROBI_003_tx 3_630166729229512303coherenceStats.mat';
path = fullfile(folder,filename);
load(path);
o1o2 = channelPairs(end-1).coherence;
end
names = {channelPairs.label};
czc4 = channelPairs(find(strcmp(names, 'Cz-C4'))).coherence;
windowSize = 1000;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
x = (1:size(czc4,1)) ./ 128;
for i = 1:size(o1o2,2)
    d = filtfilt(b, a, czc4(:,i));
    filtczc4(:,i) = d;
end
close all;
figure;
plot(x, filtczc4);
freqLabels = [{'delta'},{'theta'},{'alpha'},{'beta'},{'hibeta'}];
legend(freqLabels);
zoom xon;
pan xon;