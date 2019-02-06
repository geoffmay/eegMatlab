

%if false, adds channels one at a time, and includes all coherence pairs
individualCoherence = false;
outFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\boldEeg\neuralNets\sweep hidden layer size';
fourierFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\boldEeg\fourier\GeoffTestEEG2-edf';
convolvedFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\boldEeg\convolved';
eegFolder = 'C:\Vision\Raw Files\Geoff EEG test\history';

% %pick up where we left off
% outFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\MrEegLink';
% existingFiles = dir(outFolder);
% existingFiles([existingFiles.isdir]) = [];
% existingFiles(cellfun(@length, strfind({existingFiles.name}, 'networkDisribution')) == 0) = [];
% if(length(existingFiles) > 0)
%     latestFile = existingFiles(find([existingFiles.datenum] == max([existingFiles.datenum])));
%     [sortDates, sortInd] = sort([existingFiles.datenum]);
%     sortFiles = existingFiles(sortInd);
%     for(i = 1:length(sortFiles))
%         temp = load(fullfile(outFolder, sortFiles(i).name));
%         tempScore.training = mean([temp.best.Regressed.training]);
%         tempScore.test = mean([temp.best.Regressed.test]);
%         tempScore.validation = mean([temp.best.Regressed.validation]);
%         scores(i) = tempScore;
% %         scores(i) = temp.best.Score;
%     end
%     figure;
%     plot([[scores.training]', [scores.test]', [scores.validation]']);
%     legend({'training', 'test', 'validation'});
%     savedData = load(fullfile(outFolder, latestFile.name));
% end

%load the data
if(~exist('trainingData', 'var'))
    edfFilename = 'GeoffTestEEG2-edf.edf';
    trainingData = getConvolvedEegNeuralNetworkData(edfFilename);
    locations = fileLocations;
    
    %debug
    if(false)
        d = diff(eeg.data(1,:));
        plot(eeg.times(1:end-1), d)
        pan xon;
        zoom xon;
    end
    %end debug
end


if(individualCoherence)
    % keep = cellfun(@length, strfind(trainingData.eegColumnLabels, 'ower')) > 0;
    oldKeep = false(size(trainingData.eegColumnLabels));
    
    outputData = trainingData.mriOutput';
    
    %add one eeg parameter at a time, the one with the highest test R value
    loopCounter = 1;
    while(loopCounter < length(oldKeep))
        maxR = 0;
        maxI = -1;
        clear maxTrial;
        for i = 1:length(oldKeep)
            %         fprintf('%s: loop %d, training %d of %d (maxI = %d, maxR = %f)\n', char(datetime), loopCounter, i, length(oldKeep), maxI, maxR);
            if(~oldKeep(i))
                keep = oldKeep;
                keep(i) = true;
                inputData = trainingData.eegInput(:, keep)';
                trial.r = neuralNetworkAdvancedGeneration(inputData, outputData);
                score = trial.r.test + trial.r.validation + 0.1 * trial.r.training;
                if(score > maxR)
                    trial.keep = find(keep);
                    trial.newIndex = i;
                    maxR = score;
                    maxI = i;
                    maxTrial = trial;
                end
            end
        end
        trials(loopCounter) = maxTrial;
        oldKeep = false(size(trainingData.eegColumnLabels));
        oldKeep(maxTrial.keep) = true;
        outFile = fullfile(outFolder, sprintf('sparseEegMriNeuralNet%d.mat', loopCounter));
        save(outFile, 'trials', '-v7.3');
        fprintf('%s: loop %d, maxR = %f, measure = %s\n', char(datetime), loopCounter, maxR, trainingData.eegColumnLabels{maxTrial.newIndex});
        loopCounter = loopCounter + 1;
    end
    
    inputData = trainingData.eegInput(:, keep)';
    outputData = trainingData.mriOutput';
    regressed = neuralNetworkAdvancedGeneration(inputData, outputData);
    
    rs = [trials.r];
    plot([[rs.training]', [rs.test]', [rs.validation]']);
    legend({'training', 'test', 'validation'});
else % individualChannels == false
    distributionSize = 21;
    outputData = trainingData.mriOutput';
    inputData = trainingData.eegInput';
    %
    %     %get unique channel names
    %     channelItemIndex = 1;
    %     measureItemIndex = 2;
    %     allChans = cell(0);
    %     for i = 1:length(trainingData.eegColumnLabels)
    %         label = trainingData.eegColumnLabels{i};
    %         items = strsplit(label, ' ');
    %         if(i == 1)
    %             if(strcmp(items{1}, 'coherence') || (length(strfind(items{1}, 'ower')) > 0))
    %                 channelItemIndex = 2;
    %                 measureItemIndex = 1;
    %             elseif(strcmp(items{2}, 'coherence') || (length(strfind(items{2}, 'ower')) > 0))
    %                 channelItemIndex = 1;
    %                 measureItemIndex = 2;
    %             else
    %                 error('unhandled eeg label format: %s', label);
    %             end
    %         end
    %         chans = strsplit(items{channelItemIndex}, '-');
    %         allChans = unique([allChans, chans]);
    %     end
    %
    %     oldKeep = [];
    %     if(exist('savedData', 'var'))
    %         oldKeep = savedData.bestChan;
    %     end
    %     powerIndex = cellfun(@length, strfind(trainingData.eegColumnLabels, 'power')) > 0;
    %     coherenceIndex = cellfun(@length, strfind(trainingData.eegColumnLabels, 'coherence')) > 0;
    %
    %     for h = (length(oldKeep)+1):length(allChans)
    %         bestScore = 0;
    %         clear best loops
    %         loopCounter = 1;
    %         for i = 1:length(allChans)
    %             keepChan = oldKeep;
    %             if(~any(keepChan == i))
    %
    %                 %get input eeg data that includes exclusively channels of
    %                 %interest
    %                 keepChan(end + 1) = i;
    %                 keepMeasure = false(size(trainingData.eegColumnLabels));
    %                 for j = 1:length(keepChan)
    %                     jLabel = allChans{keepChan(j)};
    %                     hasJ = cellfun(@length, strfind(trainingData.eegColumnLabels, jLabel)) > 0;
    %                     keepMeasure(powerIndex & hasJ) = true;
    %                     for k = (j+1):length(keepChan)
    %                         kLabel = allChans{keepChan(k)};
    %                         hasK = cellfun(@length, strfind(trainingData.eegColumnLabels, kLabel)) > 0;
    %                         keepMeasure(coherenceIndex & hasJ & hasK) = true;
    %                     end
    %                 end
    %                 inputData = trainingData.eegInput(:, keepMeasure)';
    
    maxHiddenNeurons = 100;
    for h = 1:maxHiddenNeurons
        outputFilename = fullfile(outFolder, sprintf('network %d hidden neurons.mat', h));
        
        if(~exist(outputFilename, 'file'))
            
            %generate a distribution of trials
            clear regressed;
            for j = 1:distributionSize
                fprintf('%s: hidden size %d fo %d, iteration %d of %d', char(datetime), h, maxHiddenNeurons, j, distributionSize);
                tic;
                [regressed(j), nets{j}] = neuralNetworkAdvancedGeneration(inputData, outputData, h);
                times(j) = toc;
                quickScores(j) = regressed(j).training * .1 + regressed(j).test * .45 + regressed(j).validation * .45;
                fprintf(' score: %0.3f\n', quickScores(j));
            end
            medScore = median(quickScores);
            net = nets{find(quickScores == medScore)};            
            %         if(medScore > bestScore)
            %             bestChan = keepChan;
            %             bestScore = medScore;
            %             best.Score = bestScore;
            %             best.Inputs = allChans(keepChan);
            %             best.Measures = trainingData.eegColumnLabels(keepMeasure);
            %             best.Regressed = regressed;
            %         end
            loop.net = net;
            loop.score = medScore;
            loop.hiddenLayerCount = h;
            loop.secondsToComplete = times(find(quickScores == medScore));
            save(outputFilename, 'loop', '-v7.3');
            clear loop
        end
    end %h loop
end %individual coherence


if(false)
    %make a topoplot of the electrodes in order
    %topoplot(
end





