
%things to set
folder = 'E:\pca\';
folder = '/media/Seagate Backup Plus Drive/pca/';
filenameFilter = '.coherence.mat';
maxRam = 8 * 1024 * 1024 * 1024;
maxChannel = 31;

%get pertinent filenames
files =dir(folder);
files([files.isdir]) = [];
files(cellfun(@length, strfind({files.name}, filenameFilter)) == 0) = [];

%compute how much of each file we will sample
freqLabels = {'delta', 'theta', 'alpha', 'beta', 'hiBeta'};
totalFileSize = sum([files.bytes]);
sampleSkip = ceil(totalFileSize / maxRam);
columnCount = maxChannel*(maxChannel-1)/2*length(freqLabels);
if(sampleSkip > 1)
    rowCount = floor(maxRam / 8 / columnCount);
else
    rowCount = floor(totalFileSize / 8 / columnCount);
end
allData = NaN(rowCount, columnCount);
sourceIndex = NaN(rowCount, 1);

writeIndex = 1;



for fileCounter = 1:length(files)
    %load the file
    fprintf('\n%s: file %d of %d', char(datetime), fileCounter, length(files));
    fullPath = fullfile(folder, files(fileCounter).name);
    subject.name = files(fileCounter).name;
    a = load(fullPath);
    
    %compute stats for the individual subject
    for i = 1:length(a.coh.channelPairs)
        coh.label = a.coh.channelPairs(i).label;
        coh.meanCoherence = mean(a.coh.channelPairs(i).coherence, 1);
        coh.stdCoherence = std(a.coh.channelPairs(i).coherence, 1);
        subject.coherence(i) = coh;
        clear coh;
    end
    for i = 1:length(a.coh.channels)
        coh.label = a.coh.channels(i).label;
        coh.meanAbsPower = mean(a.coh.channels(i).absolutePower, 1);
        coh.stdCoherence = std(a.coh.channels(i).absolutePower, 1);
        subject.power(i) = coh;
        clear coh;
    end    
    fileCursor = 1;
    readIndex = 1;
    sampleCount = length(a.coh.x);
    
    
    %drop channels whose index is greater than maxChannel
    pairWriteCounter = 1;
    pairReadCounter = 1;
    chanReadCounter = 1;
    chanWriteCounter = 1;
    includeChannels = logical(zeros(1, length(a.coh.channels)));
    includeChannelPairs = logical(zeros(1, length(a.coh.channelPairs)));
    for i = 1:length(a.coh.channels)
        for j = (i+1):length(a.coh.channels)
            if((i <= maxChannel) & (j <= maxChannel))
                for k = 1:length(freqLabels)
                    pairLabels{pairWriteCounter} = sprintf('%s %s', a.coh.channelPairs(pairReadCounter).label, freqLabels{k});
                    pairWriteCounter = pairWriteCounter + 1;
                end
                includeChannelPairs(pairReadCounter) = 1;
            end
            pairReadCounter = pairReadCounter + 1;
        end
        if(i < maxChannel)
            for k = 1:length(freqLabels)
                chanLabels{chanWriteCounter} = sprintf('%s %s', a.coh.channels(i).label, freqLabels{k});
                chanWriteCounter = chanWriteCounter + 1;
            end
            includeChannels(chanReadCounter) = 1;
        end
        chanReadCounter = chanReadCounter + 1;
    end

    %read the data
    thisFileStart = writeIndex;
    includePairIndexes = find(includeChannelPairs);    
    while(fileCursor < sampleCount)
        counter = 1;
        for i = 1:length(includePairIndexes)
            allData(writeIndex, counter:counter + length(freqLabels)-1) = ...
                a.coh.channelPairs(includePairIndexes(i)).coherence(fileCursor, :);
            counter = counter + length(freqLabels);
        end
        sourceIndex(writeIndex) = fileCounter;
        writeIndex = writeIndex + 1;
        fileCursor = fileCursor + sampleSkip;
    end
    
    %do pca on individual
    [p.COEFF, p.SCORE, p.LATENT, p.TSQUARED, p.EXPLAINED, p.MU] = pca(allData(thisFileStart:writeIndex-1,:));
    subject.pca = p;
    
    for i = 1:size(allData,2)
        thisMean = mean(allData(thisFileStart:writeIndex-1,i));
        allData(thisFileStart:writeIndex-1,i) = allData(thisFileStart:writeIndex-1,i) - thisMean;
    end
    
    subjects(fileCounter) = subject;
    clear subject;
end
clear a;
measureLabels = pairLabels;

sampleSum = 0;
sampleSmallSum = 0;
for i = 1:length(subjects)
    sampleCount(i) = size(subjects(i).pca.SCORE, 1);
    sampleSum = sampleSum + sampleCount(i);
    smallSampleCount(i) = floor(sampleCount(i) / sampleSkip);
    sampleSmallSum = sampleSmallSum + smallSampleCount(i);
end


fprintf('\n%s saving...', char(datetime));
save(fullfile(folder, 'allWahbehMinusMeanPtsd.mat'), 'allData', 'subjects', 'measureLabels', '-v7.3');

fprintf('\n%s computing pca...', char(datetime));
[p.COEFF, p.SCORE, p.LATENT, p.TSQUARED, p.EXPLAINED, p.MU] = pca(allData);


fprintf('\n%s saving full pca...', char(datetime));
save(fullfile(folder, 'allWahbehPtsdMinusMeanPcaVerbose.mat'), 'p', 'subjects', 'measureLabels', '-v7.3');

p.SCORE = [];
p.TSQUARED = [];
fprintf('\n%s saving abbreviated pca...', char(datetime));
save(fullfile(folder, 'allWahbehPtsdMinusMeanPca.mat'), 'p', 'measureLabels', '-v7.3');

fprintf('\nfinished\n');

%compare all first components
coeff = NaN(length(subjects(1).pca.COEFF(:,1)), length(subjects));
for subjectNumber = 1:length(subjects)
    coeff(:,subjectNumber) = subjects(subjectNumber).pca.COEFF(:,1);
    explained(subjectNumber) = subjects(subjectNumber).pca.EXPLAINED(1);
    latent(subjectNumber) = subjects(subjectNumber).pca.LATENT(1);
end

Zmeasure = linkage(coeff');

Z = linkage(coeff);
figure;
dend = dendrogram(Z);

%compare means
means = NaN(length(subjects(1).pca.COEFF(:,1)), length(subjects));
stds = NaN(length(subjects(1).pca.COEFF(:,1)), length(subjects));
for subjectNumber = 1:length(subjects)
    subject = subjects(subjectNumber);
    counter = 1;
    for pairNumber = 1:length(subject.coherence)
        coh = subject.coherence(pairNumber);
        means(counter:counter+4, subjectNumber) = coh.meanCoherence;
        stds(counter:counter+4, subjectNumber) = coh.stdCoherence;
        counter = counter + 5;
    end
end

subjectSummaries.means = means;
subjectSummaries.stds = stds;
subjectSummaries.coeff1s = coeff;
save('wahbehSubjectSummaries.mat', 'subjectSummaries');


Z = linkage(means);
figure;
dend = dendrogram(Z);

Z = linkage(stds);
figure;
dend = dendrogram(Z);


Z = linkage(coeff);
figure;
dend = dendrogram(Z);

dotProducts = NaN(1, length(subjects));


fprintf('\nstarting tsne... (%s)', char(datetime));
tsneResult = tsne(coeff')
figure;
scatter(tsneResult(:,1), tsneResult(:,2))
fprintf('\nfinshed tsne... (%s)', char(datetime));

Y = tsne(coeff);

plotCoherencePca(p.COEFF(:,1), pairLabels)