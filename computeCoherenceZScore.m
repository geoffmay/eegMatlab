frequencyOfChoice = 4;
filter = '_013';
eyes = 'closed';

if(false)
    baseFilename = '/home/data/EEG/processed/Robi/coherenceReref/ROBI_003_baseline eyes open_630158995692243270coherenceStats.mat';
    exitFilename = '/home/data/EEG/processed/Robi/coherenceReref/ROBI_003_outcome eyes open_630230539149591228coherenceStats.mat';
    coherenceStats{1} = load(baseFilename);
    coherenceStats{2} = load(exitFilename);
    fileLabels = [{'pre'},{'post'}];
else
    folder = '/home/data/EEG/processed/Robi/coherenceReref/';
    files = dir(folder);
    files([files.isdir]) = [];
    hit = find(cellfun(@length, strfind({files.name}, filter)));
    
    hitFilenames = files(hit);
    baseIndex = find(cellfun(@length, strfind({hitFilenames.name}, sprintf('baseline eyes %s', eyes))));
    exitIndex = find(cellfun(@length, strfind({hitFilenames.name}, sprintf('outcome eyes %s', eyes))));
    baseFilename = fullfile(folder, files(hit(baseIndex)).name);
    exitFilename = fullfile(folder, files(hit(exitIndex)).name);
    coherenceStats{1} = load(baseFilename);
    coherenceStats{2} = load(exitFilename);
    fileLabels = [{'pre'},{'post'}];
        
    redo = false;
    if(exist('processedFilter', 'var'))
        if(strcmp(processedFilter, filter))
            redo = false;
        end
    end
    if(redo)
        clear coherenceStats fileLabels;
        counter = 1;
        for i = 1:length(hit)
            if(files(hit(i)).name(1) ~= '.')
                thisName = fullfile(folder,files(hit(i)).name);
                data = load(thisName);
                if(~isfield(data, 'hasNoise'))
                    baseFolder = '/home/data/EEG/data/ROBI/';
                    baseFile = sprintf('%s%s', baseFolder, ...
                        data.stat.filename(length('/home/data/EEG/data/ROBI')+1:end));
                    if(exist(baseFile, 'file'))
                        [hasNoise, noiseRatio] = checkForMainsNoise({baseFile});
                        data.hasNoise = hasNoise;
                        data.noiseRatio = noiseRatio;
                        stat = data.stat;
                        save(thisName, 'stat', 'hasNoise', 'noiseRatio');
                        if(~hasNoise)
                            coherenceStats{counter} = data;
                            fileLabels{counter} = files(hit(i)).name;
                            counter = counter + 1;
                        end
                    end
                else
                    if(~data.hasNoise)
                        coherenceStats{counter} = data;
                        fileLabels{counter} = files(hit(i)).name;
                        counter = counter + 1;
                    end
                end
            end
        end
        processedFilter = filter;
    end
end




agesToDetermineChannelOrder = [17 21 21 21 21 21 21 21 21 21];
freqLabels = [{'delta'},{'theta'},{'alpha'},{'beta'},{'hibeta'},{'alpha1'},{'alpha2'},{'beta1'},{'beta2'},{'beta3'}];
frequencyLookupLabels = [{'Delta 1-4Hz'}, {'Theta 4-8Hz'}, {'Alpha 8-12Hz'}, {'Beta 12-25Hz'}, {'HiBeta 25-30Hz'},{'Alpha1 8-10Hz'},{'Alpha2 10-12Hz'},{'Beta1 12-15Hz'},{'Beta2 15-18Hz'},{'Beta3 18-25Hz'}];

%load the data
filename= '/home/data/EEG/processed/Robi/zScore.mat'
load(filename);

%determine rows of interest
eyeOpen = strcmp(data{:,'Subject_Eyes'}, eyes);
absPower = strcmp(data{:,'Analytic_Type'}, 'absolute power mean');
coherence = strcmp(data{:,'Analytic_Type'}, 'coherence mean');
coherenceStd = strcmp(data{:,'Analytic_Type'}, 'coherence standard deviation');
freqRows = strcmp(data{:,'Spectral_Bin_Label'}, frequencyLookupLabels{frequencyOfChoice});

%get channel/pair labels
refRows = find(eyeOpen & absPower & freqRows);
meanRows = find(eyeOpen & coherence & freqRows);
stdRows = find(eyeOpen & coherenceStd & freqRows);

channelPairLabels = data{meanRows, 'Channel_Pair_Label'};
channelLabels = data{refRows, 'Channel_Pair_Label'};
searchLabels = channelLabels;
for i = 1:length(searchLabels)
    searchLabels = strrep(searchLabels, 'T3', 'T7');
    searchLabels = strrep(searchLabels, 'T4', 'T8');
    searchLabels = strrep(searchLabels, 'T5', 'P7');
    searchLabels = strrep(searchLabels, 'T6', 'P8');
end
channelIndexes = NaN(1, length(searchLabels));
for i = 1:length(channelIndexes)
    channelIndexes(i) = find(strcmp(coherenceStats{1}.stat.channelLabels, searchLabels{i}));
end

ageIndex = 20; %age 30 or so
mat = zeros(length(channelIndexes), length(channelIndexes), length(coherenceStats));

%build matrix
for fileIndex = 1:length(coherenceStats)
    counter = 1;
    for chan1 = 1:length(channelLabels)
        for chan2 = (chan1 + 1):length(channelLabels)
            thisMean = data{meanRows(counter), ageIndex};
            thisStd = data{stdRows(counter), ageIndex};
            sourceChan1 = channelIndexes(chan1);
            sourceChan2 = channelIndexes(chan2);
            if(sourceChan2 < sourceChan1)
                temp = sourceChan1;
                sourceChan1 = sourceChan2;
                sourceChan2 = temp;
            end
            thisCoh = coherenceStats{fileIndex}.stat.meanCoherences(sourceChan1, sourceChan2, frequencyOfChoice);
            %        thisCoh = baseline{fileIndex}.stat.stdDevCoherences(sourceChan1, sourceChan2, frequencyOfChoice);
            thisCoh = thisCoh * 100;
            %debug
            %end debug
            zScore = (thisCoh - thisMean) / thisStd;
            mat(chan1,chan2, fileIndex) = zScore;
            mat(chan2,chan1, fileIndex) = zScore;
            counter = counter + 1;
        end
    end
end

%[sortMat, masterPerms] = clusterMatrix(mat(:,:,end));

%order from beta, eldest group
betaOptimalPerm = [8 19 7 15 9 10 16 3 17 4 5 6 18 1 11 2 12 13 14];
%manually reordered, optimized delta across several ages
deltaOptimalPerm = [19 8 7 15 9 10 16 6 5 18 4 17 3 1 2 12 14 13 11 ];
%manually reordered, attempted compromise between beta and delta
compromisePerm = [ 8 19 7 15 9 10 16 5 18 6 4 17 3 1 2 12 14 13 11 ];

masterPerms = betaOptimalPerm;

%[sortMat, masterPerms] = clusterMatrix(mat(:,:,agesToDetermineChannelOrder(frequencyOfChoice)));

permLabels = channelLabels(masterPerms);
globalMin = min(min(min(mat)));
globalMax = max(max(max(mat)));

%plot everything
close all;
for fileCounter = 1:length(coherenceStats)
    fig = figure;
    colormap(fig, parula);
    [sortMat] = clusterMatrix(mat(:,:,fileCounter),masterPerms);
    imagesc(sortMat, [globalMin, globalMax]);
    colorbar;
    fileLabel = fileLabels{fileCounter};
    %  ageLabel = data.Properties.VariableNames{ageCounter+6};
    title(sprintf('%s %s', strrep(fileLabel, '_', ' '), freqLabels{frequencyOfChoice}));
    % sortLabels = channelLabels(perms(ageCounter,:));
    set(gca,'xtick',1:length(permLabels));
    set(gca,'ytick',1:length(permLabels));
    set(gca,'xticklabels',permLabels);
    set(gca,'yticklabels',permLabels);
    
    counter = 1;
    for i = 1:size(mat, 1)
        for j = i+1:size(mat, 2)
            pairLabels{counter} = sprintf('%s-%s', channelLabels{i}, channelLabels{j});
            pairValues(counter,:) = mat(i, j, :);
            counter = counter + 1;
        end
    end
    cMin = min(min(min(mat)));
    cMax = max(max(max(mat)));
    extreme = max(abs(cMin), abs(cMax));
    plotChannelPairs(pairLabels, pairValues(:,fileCounter), [-extreme, extreme]);
    title(sprintf('%s %s %s', strrep(fileLabel, '_', ' '), freqLabels{frequencyOfChoice}, filter));
end


tilefigs;
