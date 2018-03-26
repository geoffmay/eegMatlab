doBigPlot = 1;

folder = '/home/data/EEG/processed/Robi/coherencePca';
labelFilename = 'labels.mat';
load(fullfile(folder,labelFilename));
labels = allLabels;

files = dir(folder);
files([files.isdir]) = [];
files(find(strcmp({files.name}, labelFilename)))= [];

if(~exist('eeglab','file'))
    addpath('/home/data/EEG/scripts/eeglab13_4_4b')
    eeglab;
end
close all;

clear grandSum;

for i = 1:length(files)
    path = fullfile(folder,files(i).name);
    load(path);
    coefficient = 1;
    prefix = sprintf('%s comp %d (%.1f%%)',files(i).name, coefficient, cohPca.EXPLAINED(coefficient));
    prefix = strrep(prefix, '_', ' ');
    if(doBigPlot)
        plotCoherencePca(cohPca, labels, coefficient, prefix);
    end
    
    score = cohPca.SCORE;
    sampleRate = 128;
    windowSize = 100;
    interval = 1 * sampleRate;
    
    startFrame = 1;
    endFrame = startFrame + sampleRate*windowSize - 1;
    fftSum = zeros(1,floor(endFrame/50));
    counter = 0;
    x = ((1:length(fftSum)) ./ windowSize);
    
    while(endFrame <= size(score,1))
        window = score(startFrame:endFrame);
        scoreFft = fft(window);
        modFft = abs(scoreFft);
        modFft = modFft(1:length(fftSum));
        modFft = modFft .* x;
        fftSum = fftSum + modFft;
        counter = counter + 1;
        startFrame = startFrame + interval;
        endFrame = startFrame + sampleRate*windowSize - 1;        
    end
    fftAvg = fftSum ./ counter;
    figure;
    plot(x, fftAvg);
    file = strrep(files(i).name, '_', ' ');
    
    title(sprintf('pca fourier %s', file));
    
    if(~exist('grandSum','var'))
        grandSum = fftAvg;
    else
        grandSum = grandSum + fftAvg;
    end
       

end
grandAvg = grandSum ./ length(files);

figure;
plot(x, grandAvg);
file = 'grand average';

title(sprintf('pca fourier %s', file));

tilefigs;