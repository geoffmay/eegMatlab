%flags
doPlot = 0;
close all;

%filenames
inputFolder1 = '/media/eegDrive';
outputFolder = '/home/data/EEG/processed/Robi/zCorrelation';
inputFiles = dir(inputFolder1);

%filter
inputFiles([inputFiles.isdir]) = [];
hasTx = cellfun(@length, strfind({inputFiles.name}, 'tx'));
hasImp = cellfun(@length, strfind({inputFiles.name}, 'impedance'));
remove = find(~hasTx | hasImp);
inputFiles(remove) = [];

% clear coherenceSum;
% clear coherenceCounter;

for fileCounter = 1:length(inputFiles)
    fprintf('\n%s (%d of %d)', char(datetime), fileCounter, length(inputFiles));
    coherenceFile = inputFiles(fileCounter).name;
    shortName = coherenceFile(1:strfind(coherenceFile, '_63')-1);
    outputFilename = fullfile(outputFolder, sprintf('%s.mat',shortName));
    if(~exist(outputFilename, 'file'))
        try
            %load coherence data
            fprintf('loading...');
            cohData = load(fullfile(inputFolder1, coherenceFile));
            if(isfield(cohData, 'channelPairs'))
                sampleCount = size(cohData.channelPairs(1).coherence,1);
                cohTimes = .5:(1/128):(.5 + (sampleCount-1) / 128);
                rawDataFile = cohData.filename;
                rawDataFile = strrep(rawDataFile, '.eegData', 'Events.txt');
                
                %read the (raw) file
                fprintf('loading...');
                text = textread(rawDataFile, '%s');
                numColumns = 13;
                numRows = length(text) / numColumns;
                text1 = text;
                if(numRows ~= floor(numRows))
                    numRows = floor(numRows);
                    text1 = text1(1:(numRows * numColumns));
                end
                text2 = reshape(text1, [numColumns, numRows])';
                
                %parse the text
                numeric = zeros(size(text2, 2));
                numbers = NaN(size(text2));
                fprintf('parsing......');
                for pairCounter = 1:size(text2, 2)
                    fprintf('%d', pairCounter);
                    try
                        a = str2double(text2(1, pairCounter));
                        numeric(pairCounter) = 1;
                    catch
                    end
                    if(numeric(pairCounter))
                        for j = 1:size(text2, 1)
                            debug = false;
                            if(debug)
                                quickText = text2{j,pairCounter};
                                quickNum1 = parseDecimal(quickText);
                                quickNum2 = str2num(quickText);
                                if(quickNum1 ~= quickNum2)
                                    error(sprintf('inconsistent parsing: \n%s\n%.100f\n%.100f', quickText, quickNum1, quickNum2));
                                end
                            end
                            numbers(j,pairCounter) = parseDecimal(text2{j,pairCounter});
                        end
                    end
                end
                elapsed = toc;
                rawTimes = numbers(:,1);
                smoothZScores = numbers(:, 3);
                momentaryZScores = numbers(:, 5);
                squareSizes = numbers(:, 7);
                positiveFeedback = numbers(:, 9);
                powerZScore = numbers(:, 11);
                coherenceZScore = numbers(:, 13);
                
                
                timeCourse.rawTimes = numbers(:,1);
                timeCourse.smoothZScores = numbers(:, 3);
                timeCourse.momentaryZScores = numbers(:, 5);
                timeCourse.squareSizes = numbers(:, 7);
                timeCourse.positiveFeedback = numbers(:, 9);
                timeCourse.powerZScore = numbers(:, 11);
                timeCourse.coherenceZScore = numbers(:, 13);
                
                
                %unwrap times that have overflowed
                fprintf('processing...');
                processorFrequency = 2530957;
                processorReciprocal = 1 / processorFrequency;
                wrapValue = power(2,32) / processorFrequency;
                a = diff(numbers(:,1));
                wrapPoints = find(abs(a) > 800);
                for pairCounter = 1:length(wrapPoints)
                    rawTimes(wrapPoints(pairCounter)+1:end) = rawTimes(wrapPoints(pairCounter)+1:end) + wrapValue;
                end
                
                %remove nans
                remove = isnan(rawTimes);
                a = diff(rawTimes);
                decreasing = find(a < 0);
                if(length(decreasing) > 0)
                    remove(decreasing(1):end) = 1;
                end
                rawTimes(remove) = [];
                smoothZScores(remove) = [];
                momentaryZScores(remove) = [];
                squareSizes(remove) = [];
                positiveFeedback(remove) = [];
                powerZScore(remove) = [];
                coherenceZScore(remove) = [];
                %interpolate values
                interpTotalZ = interp1(rawTimes, momentaryZScores, cohTimes);
                interpSmoothZ = interp1(rawTimes, smoothZScores, cohTimes);
                keep = ~isnan(interpTotalZ);
                
                %determine whether references have dropped
                correlation.badRefIndex = checkForBadReference(cohData.channelPairs(31).coherence);
                keep2 = keep;
                if(correlation.badRefIndex < length(keep2))
                    keep2(correlation.badRefIndex:end) = 0;
                end
                
                %correlate
                fprintf('xxx');
                for pairCounter = 1:length(cohData.channelPairs)
                    fprintf('\b\b\b%03d', pairCounter)
                    cohPlot = cohData.channelPairs(pairCounter).coherence;
                    for freqCounter = 1:size(cohPlot,2)
                        coh = cohPlot(:,freqCounter);
                        [correlation.rho(pairCounter, freqCounter), correlation.p(pairCounter, freqCounter)] = corr(coh(keep), interpTotalZ(keep)');
                        [correlation.smoothRho(pairCounter, freqCounter), correlation.smoothP(pairCounter, freqCounter)] = corr(coh(keep), interpSmoothZ(keep)');
                        if(length(find(keep2)) > (2048 * 10))
                            [correlation.goodRho(pairCounter, freqCounter), correlation.goodP(pairCounter, freqCounter)] = corr(coh(keep2), interpTotalZ(keep2)');
                            [correlation.goodSmoothRho(pairCounter, freqCounter), correlation.goodSmoothP(pairCounter, freqCounter)] = corr(coh(keep2), interpSmoothZ(keep2)');
                        else
                            [correlation.goodRho(pairCounter, freqCounter), correlation.goodP(pairCounter, freqCounter)] = NaN;
                            [correlation.goodSmoothRho(pairCounter, freqCounter), correlation.goodSmoothP(pairCounter, freqCounter)] = NaN;
                        end
                        correlation.meanValues(pairCounter, freqCounter) = mean(coh(keep));
                        correlation.stdValues(pairCounter, freqCounter) = std(coh(keep));
                        correlation.meanGoodValues(pairCounter, freqCounter) = mean(coh(keep2));
                        correlation.stdGoodValues(pairCounter, freqCounter) = std(coh(keep2));
                    end
                end
                correlation.meanZ = interpTotalZ(keep);
                correlation.stdZ = interpTotalZ(keep);
                correlation.goodMeanZ = interpTotalZ(keep2);
                correlation.goodStdZ = interpTotalZ(keep2);
                
                %plot
                correlation.channelPairLabels = {cohData.channelPairs.label};
                correlation.frequencyLabels = [{'delta'},{'theta'},{'alpha'},{'beta'},{'hibeta'}];
                if(doPlot)
                    figure;
                    imagesc(correlation.smoothRho);
                    colorbar;
                    myTitle = sprintf('smooth %s', shortName);
                    myTitle = strrep(myTitle, '_', ' ');
                    title(myTitle);
                    threshold = 95;
                    supraThreshold = prctile( correlation.goodSmoothRho, threshold)
                    subThreshold = prctile( correlation.goodSmoothRho, 100 - threshold)
                    maxFreq = find(supraThreshold == max(supraThreshold));
                    minFreq = find(subThreshold == min(subThreshold));
                    maxPair = correlation.goodSmoothRho(:, maxFreq) > supraThreshold(maxFreq);
                    minPair = correlation.goodSmoothRho(:, maxFreq) > supraThreshold(maxFreq);
                    plotChannelPairs(correlation.channelPairLabels(maxPair));
                    correlation.maxCorrelatedFrequency = correlation.frequencyLabels{maxFreq};
                    correlation.maxCorrelatedChannelPairs = correlation.channelPairLabels(maxPair);
                    myTitle = sprintf('smooth %s %s', shortName, correlation.frequencyLabels{maxFreq});
                    myTitle = strrep(myTitle, '_', ' ');
                    title(myTitle);
                    drawnow;
                end
                
                %save
                correlation.filename = rawDataFile;
                save(outputFilename, 'correlation', 'timeCourse');
            end
        catch err
            save(outputFilename, 'err');
            %             rethrow(err);
        end
    else
%         if(~exist('coherenceSum', 'var'))
%             coherenceSum = correlation.rho .* correlation.badRefIndex;
%             coherenceCount = correlation.badRefIndex;
%         else
%             coherenceSum = coherenceSum + correlation.rho .* correlation.badRefIndex;
%             coherenceCount = coherenceCount + correlation.badRefIndex;
%         end
        if(doPlot)
            if(false)
                data{fileCounter} = load(outputFilename);
                plotChannelPairs(data{fileCounter}.correlation.maxCorrelatedChannelPairs);
                [folder file ext] = fileparts(data{fileCounter}.correlation.filename);
                file = folder(strfind(folder, '/ROBI/') + length('/ROBI/'):end);
                myTitle = sprintf('%s %s', file, data{fileCounter}.correlation.maxCorrelatedFrequency);
                myTitle = strrep(myTitle, '_', ' ');
                title(myTitle);
                drawnow;
            else
                %plot saved output
                if(strfind(outputFilename, 'ROBI_003'))
                    data{fileCounter} = load(outputFilename);
                    if(isfield(data{fileCounter}, 'correlation'))
                        correlation = data{fileCounter}.correlation;
                        figure;
                        imagesc(correlation.smoothRho);
                        colorbar;
                        myTitle = sprintf('smooth %s', shortName);
                        myTitle = strrep(myTitle, '_', ' ');
                        title(myTitle);
                        if(false)
                            threshold = 95;
                            supraThreshold = prctile( correlation.goodSmoothRho, threshold)
                            subThreshold = prctile( correlation.goodSmoothRho, 100 - threshold)
                            maxFreq = find(supraThreshold == max(supraThreshold));
                            minFreq = find(subThreshold == min(subThreshold));
                            maxPair = correlation.goodSmoothRho(:, maxFreq) > supraThreshold(maxFreq);
                            minPair = correlation.goodSmoothRho(:, maxFreq) > supraThreshold(maxFreq);
                            plotChannelPairs(correlation.channelPairLabels(maxPair));
                            correlation.maxCorrelatedFrequency = correlation.frequencyLabels{maxFreq};
                            correlation.maxCorrelatedChannelPairs = correlation.channelPairLabels(maxPair);
                            myTitle = sprintf('smooth %s %s', shortName, correlation.frequencyLabels{maxFreq});
                            myTitle = strrep(myTitle, '_', ' ');
                            title(myTitle);
                        end
                        drawnow;
                    end
                end
            end
        end
    end
end


source = '/home/data/EEG/processed/Robi/fullZ/noPca/ROBI_003-baseline eyes open-630158995692243270.zScore';
dest = '/home/data/EEG/processed/Robi/fullZ/noPca/temp.mat';
copyfile(source, dest);
baseline = load(dest);

%todo: put this into a function, and remove mastoids automatically

avgCoherence = coherenceSum ./ coherenceCount;
figure;
imagesc(avgCoherence);
colorbar;
threshold = 95;
supraThreshold = prctile( avgCoherence, threshold)
subThreshold = prctile( avgCoherence, 100 - threshold)
maxFreq = find(supraThreshold == max(supraThreshold));
minFreq = find(subThreshold == min(subThreshold));
maxPair = avgCoherence(:, maxFreq) > supraThreshold(maxFreq);
minPair = avgCoherence(:, minFreq) < subThreshold(minFreq);
maxLabels = correlation.channelPairLabels(maxPair)';
minLabels = correlation.channelPairLabels(minPair)';
plotChannelPairs([{maxLabels}, {minLabels}]);
correlation.maxCorrelatedFrequency = correlation.frequencyLabels{maxFreq};
correlation.maxCorrelatedChannelPairs = correlation.channelPairLabels(maxPair);
myTitle = sprintf('smooth %s %s', 'average', correlation.frequencyLabels{maxFreq});
myTitle = strrep(myTitle, '_', ' ');
title(myTitle);

tilefigs;
