plotResult = 0;
onlyRest = 0;
saveResult = 1;
reduceDiskSize = 0;

inputFolder = '/media/eegDrive/';
inputFolder = '/media/eegDrive/diffRef/';
outputFolder = '/home/data/EEG/processed/Robi/coherencePca/';

inputFiles = dir(inputFolder);
inputFiles([inputFiles.isdir]) = [];

if(onlyRest)
    %filter
    outcome = cellfun(@length, strfind({inputFiles.name},'outcome'));
    baseline = cellfun(@length, strfind({inputFiles.name},'baseline'));
    remove = find(~baseline & ~outcome);
    inputFiles(remove) = [];
end

for fileNumber = 1:length(inputFiles)
    inputFile = inputFiles(fileNumber).name;
    inputPath = fullfile(inputFolder, inputFile);
    fprintf('\n%s file# %d of %d', char(datetime), fileNumber, length(inputFiles));
    truncated = inputFile(1:strfind(inputFile,'_630')-1);
    if(length(truncated) > 0)
        outputFile = sprintf('%s%s.mat', outputFolder, truncated);
    else
        outputFile = sprintf('%s%s', outputFolder, inputFile);
    end
    %if(~exist(outputFile, 'file'))
    if(true)
        try
            coherenceData = load(inputPath);
            if(isfield(coherenceData, 'channelPairs'))
                badReferenceFrame = checkForBadReference(coherenceData.channelPairs(31).coherence(:,1));
                if(badReferenceFrame < size(coherenceData.channelPairs(31).coherence, 1))
                    for i = 1:length(coherenceData.channelPairs)
                        coherenceData.channelPairs(i).coherence(badReferenceFrame:end,:) = [];
                    end
                end
                
                %determine output filename
                [folder, file, ext] = fileparts(coherenceData.filename);
                file = strrep(folder, '/home/data/EEG/data/ROBI/', '');
                file = strrep(file, '/', ' ');
                
                
                allSize = size(coherenceData.channelPairs(1).coherence);
                allSize(2) = allSize(2) * length(coherenceData.channelPairs);
                
                allDerivations = NaN(allSize);
                freqLabels = {'delta','theta','alpha','beta','hibeta'};
                
                start = 1;
                finish = start + 4;
                for i = 1:length(coherenceData.channelPairs)
                    allDerivations(:, start:finish) = coherenceData.channelPairs(i).coherence;
                    for(j = 1:length(freqLabels))
                        allLabels{start + j - 1} = sprintf('%s %s', coherenceData.channelPairs(i).label, freqLabels{j});
                    end
                    
                    start = start + 5;
                    finish = start + 4;
                    
                end
                singleLabels = antChannelLocs;
                singleLabels = singleLabels(1:33);
                for i = 1:length(coherenceData.channels)
                    allDerivations(:, start:finish) = coherenceData.channels(i).absolutePower(:,1:length(freqLabels));
                    for j = 1:length(freqLabels)
                        allLabels{start + j - 1} = sprintf('%s %s abs', singleLabels{i}, freqLabels{j});
                    end
                    start = start + 5;
                    finish = start + 4;
                    allDerivations(:, start:finish) = coherenceData.channels(i).relativePower(:,1:length(freqLabels));
                    for j = 1:length(freqLabels)
                        allLabels{start + j - 1} = sprintf('%s %s rel', singleLabels{i}, freqLabels{j});
                    end
                    start = start + 5;
                    finish = start + 4;
                end
                
                %remove mastoids
                m1 = cellfun(@length, strfind(allLabels, 'M1'));
                m2 = cellfun(@length, strfind(allLabels, 'M2'));
                cpz = cellfun(@length, strfind(allLabels, 'CPz'));
                remove = find(m1 | m2 | cpz);
                allDerivations(:,remove) = [];
                allLabels(:,remove) = [];
                
                %perform pca
                if(exist('eeglab','file'))
                    rmpath(genpath(fileparts(which('eeglab'))));
                end
                [cohPca.COEFF, cohPca.SCORE, cohPca.LATENT, cohPca.TSQUARED, cohPca.EXPLAINED, cohPca.MU] = pca(allDerivations);
                close all;
                
                %plot first coefficient
                if(plotResult)
                    plotCoherencePca(cohPca, allLabels, 1);
                end
                
                
                %get channel lables only (no frequency)
                onlyLabels = allLabels;
                for i = 1:length(onlyLabels)
                    label = onlyLabels{i};
                    space = strfind(label, ' ');
                    if(length(space) > 0)
                        label(space(1):end) = [];
                    end
                    onlyLabels{i} = label;
                end
                cohPca.labelChannelOnly = onlyLabels;
                cohPca.labelChannelFreq = allLabels;
                cohPca.filename = coherenceData.filename;
                
                %display image of pair coefficients
                if(false)
                    imagesc(cohPca.COEFF(:, 1:10));
                    colorbar;
                end
                
                %check top values
                doPlot = 0;
                doTopo = 0;
                maxPair = length(cohPca.labelChannelOnly);
                maxCoefficient = 20;
                
                if(doTopo)
                    [chanLabs, chanlocs] = antChannelLocs;
                else
                    [chanLabs] = antChannelLocs;
                end
                chanCoefficients = zeros(length(chanLabs), maxCoefficient);
                
                M1 = strcmp(chanLabs, 'M1');
                M2 = strcmp(chanLabs, 'M2');
                CPz = strcmp(chanLabs, 'CPz');
                remove = M1 | M2 | CPz;
                
                if(false)
                    for i = 1:maxCoefficient
                        a = cohPca.COEFF(:, i);
                        [a1, ind] = sort(a, 1, 'descend');
                        if(doPlot)
                            if(~doTopo)
                                figure;
                                plot(a1);
                            end
                        end
                        sortLabels = cohPca.labelChannelOnly(ind);
                        maxWeight = max(cohPca.COEFF(:, i));
                        alpha = cohPca.COEFF(:, i) ./ maxWeight;
                        if(false)
                            i
                            %plot sums from top channels
                            for chanPairCounter = 1:maxPair
                                chans = strsplit(sortLabels{chanPairCounter}, '-');
                                ind1 = find(strcmp(chanLabs, chans{1}));
                                ind2 = find(strcmp(chanLabs, chans{2}));
                                
                                chanCoefficients(ind1, i) = chanCoefficients(ind1, i) + cohPca.COEFF(chanPairCounter, i);
                                chanCoefficients(ind2, i) = chanCoefficients(ind2, i) + cohPca.COEFF(chanPairCounter, i);
                            end
                            if(doPlot)
                                figure;
                                topoplot(chanCoefficients(~remove,i), chanlocs(~remove), 'maplimits', [-3 3]);
                                colorbar;
                            end
                        end
                        sortLabels(1:20)
                        if(doPlot)
                            if(doTopo)
                                plotChannelPairs(sortLabels(1:100));
                                title(sprintf('component %d %s', i, cohPca.filename));
                            end
                        end
                    end
                end
                
                %reduce size to save disk space
                cohPca.COEFF(:, maxCoefficient+1:end) = [];
                cohPca.SCORE(:, maxCoefficient+1:end) = [];
                cohPca.EXPLAINED(:, maxCoefficient+1:end) = [];
                
                save(outputFile, 'cohPca', 'chanCoefficients', 'allLabels');
            elseif(isfield(coherenceData, 'timeCourse'))
                %determine output filename
                [folder, file, ext] = fileparts(coherenceData.timeCourse.basefile);
                file = strrep(folder, '/home/data/EEG/data/ROBI/', '');
                file = strrep(file, '/', ' ');
                freqLabels = {'delta','theta','alpha','beta','hibeta'};

                singleLabels = antChannelLocs;
                singleLabels = singleLabels(1:33);

                %remove mastoids
                m1 = cellfun(@length, strfind(singleLabels, 'M1'));
                m2 = cellfun(@length, strfind(singleLabels, 'M2'));
                cpz = cellfun(@length, strfind(singleLabels, 'CPz'));
                removeSingle = find(m1 | m2 | cpz);
                removeDouble = zeros(1,528);
                doubleLabels = cell(size(removeDouble));
                counter = 1;
                for i = 1:33
                    badI = length(find(removeSingle==i)) > 0;
                    for j = i+1:33
                        badJ = length(find(removeSingle==j)) > 0;
                        if(badI | badJ)
                            removeDouble(counter) = 1;
                        end
                        doubleLabels{counter} = sprintf('%s-%s', singleLabels{i}, singleLabels{j});
                        counter = counter + 1;
                    end
                end
                coherenceData.timeCourse.coherencePlot(:,:,find(removeDouble)) = [];
                coherenceData.timeCourse.powerPlot(:,:,removeSingle) = [];
                doubleLabels(find(removeDouble)) = [];
                

                %reshape data with proper labels
                cohSize = size(coherenceData.timeCourse.coherencePlot);
                powSize = size(coherenceData.timeCourse.powerPlot);
                allSize = [cohSize(1), cohSize(2) * cohSize(3) + powSize(2) * powSize(3)];
                allDerivations = NaN(allSize);
                allLabels = cell(1, allSize(2));
                freqLabels = {'delta','theta','alpha','beta','hibeta'};
                
                start = 1;
                finish = start + 4;
                for i = 1:cohSize(3)
                    allDerivations(:, start:finish) = coherenceData.timeCourse.coherencePlot(:,:,i);
                    for(j = 1:length(freqLabels))
                        allLabels{start + j - 1} = sprintf('%s %s', doubleLabels{i}, freqLabels{j});
                    end                    
                    start = start + 5;
                    finish = start + 4;                    
                end
                for i = 1:powSize(3)
                    allDerivations(:, start:finish) = coherenceData.timeCourse.coherencePlot(:,:,i);
                    for j = 1:length(freqLabels)
                        allLabels{start + j - 1} = sprintf('%s %s abs', singleLabels{i}, freqLabels{j});
                    end                    
                    start = start + 5;
                    finish = start + 4;                    
                end
                
                
                
                %                 singleLabels = antChannelLocs;
                %                 singleLabels = singleLabels(1:33);
                %                 for i = 1:length(coherenceData.channels)
                %                     allDerivations(:, start:finish) = coherenceData.channels(i).absolutePower(:,1:length(freqLabels));
                %                     for j = 1:length(freqLabels)
                %                         allLabels{start + j - 1} = sprintf('%s %s abs', singleLabels{i}, freqLabels{j});
                %                     end
                %                     start = start + 5;
                %                     finish = start + 4;
                %                     allDerivations(:, start:finish) = coherenceData.channels(i).relativePower(:,1:length(freqLabels));
                %                     for j = 1:length(freqLabels)
                %                         allLabels{start + j - 1} = sprintf('%s %s rel', singleLabels{i}, freqLabels{j});
                %                     end
                %                     start = start + 5;
                %                     finish = start + 4;
                %                 end

                
                %perform pca
                if(exist('eeglab','file'))
                    rmpath(genpath(fileparts(which('eeglab'))));
                end
                [cohPca.COEFF, cohPca.SCORE, cohPca.LATENT, cohPca.TSQUARED, cohPca.EXPLAINED, cohPca.MU] = pca(allDerivations);
                close all;
                
                %plot first coefficient
                if(plotResult)
                    plotCoherencePca(cohPca, allLabels, 1);
                end
                
                
                %get channel lables only (no frequency)
                onlyLabels = allLabels;
                for i = 1:length(onlyLabels)
                    label = onlyLabels{i};
                    space = strfind(label, ' ');
                    if(length(space) > 0)
                        label(space(1):end) = [];
                    end
                    onlyLabels{i} = label;
                end
                cohPca.labelChannelOnly = onlyLabels;
                cohPca.labelChannelFreq = allLabels;
                cohPca.filename = coherenceData.timeCourse.basefile;
                
                %display image of pair coefficients
                if(false)
                    imagesc(cohPca.COEFF(:, 1:10));
                    colorbar;
                end
                
                %check top values
                doPlot = 0;
                doTopo = 0;
                maxPair = length(cohPca.labelChannelOnly);
                maxCoefficient = 20;
                
                if(doTopo)
                    [chanLabs, chanlocs] = antChannelLocs;
                else
                    [chanLabs] = antChannelLocs;
                end
                chanCoefficients = zeros(length(chanLabs), maxCoefficient);
                
                M1 = strcmp(chanLabs, 'M1');
                M2 = strcmp(chanLabs, 'M2');
                CPz = strcmp(chanLabs, 'CPz');
                remove = M1 | M2 | CPz;
                
                
                
                %                 if(false)
                %                     for i = 1:maxCoefficient
                %                         a = cohPca.COEFF(:, i);
                %                         [a1, ind] = sort(a, 1, 'descend');
                %                         if(doPlot)
                %                             if(~doTopo)
                %                                 figure;
                %                                 plot(a1);
                %                             end
                %                         end
                %                         sortLabels = cohPca.labelChannelOnly(ind);
                %                         maxWeight = max(cohPca.COEFF(:, i));
                %                         alpha = cohPca.COEFF(:, i) ./ maxWeight;
                %                         if(false)
                %                             i
                %                             %plot sums from top channels
                %                             for chanPairCounter = 1:maxPair
                %                                 chans = strsplit(sortLabels{chanPairCounter}, '-');
                %                                 ind1 = find(strcmp(chanLabs, chans{1}));
                %                                 ind2 = find(strcmp(chanLabs, chans{2}));
                %
                %                                 chanCoefficients(ind1, i) = chanCoefficients(ind1, i) + cohPca.COEFF(chanPairCounter, i);
                %                                 chanCoefficients(ind2, i) = chanCoefficients(ind2, i) + cohPca.COEFF(chanPairCounter, i);
                %                             end
                %                             if(doPlot)
                %                                 figure;
                %                                 topoplot(chanCoefficients(~remove,i), chanlocs(~remove), 'maplimits', [-3 3]);
                %                                 colorbar;
                %                             end
                %                         end
                %                         sortLabels(1:20)
                %                         if(doPlot)
                %                             if(doTopo)
                %                                 plotChannelPairs(sortLabels(1:100));
                %                                 title(sprintf('component %d %s', i, cohPca.filename));
                %                             end
                %                         end
                %                     end
                %                 end
                
                if(reduceDiskSize)
                    %reduce size to save disk space
                    cohPca.COEFF(:, maxCoefficient+1:end) = [];
                    cohPca.SCORE(:, maxCoefficient+1:end) = [];
                    cohPca.EXPLAINED(:, maxCoefficient+1:end) = [];
                end
                
                if(saveResult)                
                    save(outputFile, 'cohPca', 'chanCoefficients', 'allLabels');
                end
            end
            
        catch err
            output.err = err;
            output.date = char(datetime);
            save(outputFile, 'output');
        end
        
    end
end

