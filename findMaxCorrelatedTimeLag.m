windowDuration = 0.5;
windowInterval = 0.25;
searchLimits = [-0.25, 0.25];
searchFrameInterval = 1;
% filename = 'C:\Vision\Raw Files\eegtest\export\eegtest_0001-edf.edf';
% markerPath = 'C:\Vision\Raw Files\eegtest\export\eegtest_0001_Pulse Artifact Correction.Markers';
% outputFolder = 'C:\Vision\Raw Files\eegtest\export\corrLag';
% filename = 'C:\Vision\Raw Files\Geoff EEG test\export\GeoffTestEEG2-edf.edf';
% markerPath = 'C:\Vision\Raw Files\Geoff EEG test\export\GeoffTestEEG2_Pulse Artifact Correction.Markers';
outputFolder = 'C:\Vision\Raw Files\Geoff EEG test\export\corrLag';

inputFolder = 'C:\Vision\Raw Files\Geoff EEG test\history\';
inputFiles = dir(inputFolder);
inputFiles = inputFiles(cellfun(@length, strfind({inputFiles.name}, '-edf.edf')) > 0);

for fileCounter = 2:length(inputFiles)
    inputPath = fullfile(inputFolder, inputFiles(fileCounter).name);
    markerPath = strrep(inputPath, '-edf.edf', '_Pulse Artifact Correction.Markers');
    filename = inputPath;
    
    outputSubFolder = fullfile(outputFolder, inputFiles(fileCounter).name);
    %     if(~exist(outputSubFolder))
    mkdir(outputSubFolder);
    
    
    eeg = loadBrainvisionEdf(filename);
    ecgChannelI = find(cellfun(@length, strfind({eeg.chanlocs.labels}, 'ECG'))>0);
    if(length(ecgChannelI) == 1)
        eeg.data(ecgChannelI, :) = [];
        eeg.nbchan = eeg.nbchan - 1;
        eeg.chanlocs(ecgChannelI) = [];
    end
    positions = loadVolumeMarkers(markerPath);
    
    windowCount= eeg.pnts / (windowInterval * eeg.srate) - 1;
    
    searchRange = ceil(searchLimits(1)*eeg.srate) : searchFrameInterval : floor(searchLimits(2)*eeg.srate);
    startIndex1 = abs(min(searchRange)) + 1;
    endIndex1 = startIndex1 + (windowDuration) * eeg.srate - 1;
    
    % lagMatrix = NaN(eeg.nbchan, eeg.nbchan, windowCount);
    % rhoMatrix = NaN(eeg.nbchan, eeg.nbchan, windowCount);
    % pMatrix = NaN(eeg.nbchan, eeg.nbchan, windowCount);
    for windowCounter = 1:windowCount
        outputFile = fullfile(outputSubFolder, sprintf('corrLag-%06d.mat', windowCounter));
        if(~exist(outputFile, 'file'))
            placeHolder = sprintf('started on %s', char(datetime));
            save(outputFile, 'placeHolder');
            fprintf('\n%s: window %d of %d', char(datetime), windowCounter, windowCount);
            startIndex1 = abs(min(searchRange)) + 1 + floor((windowCounter-1) * windowInterval * eeg.srate);
            endIndex1 = startIndex1 + (windowDuration) * eeg.srate - 1;
            lagMatrix = NaN(eeg.nbchan, eeg.nbchan);
            rhoMatrix = NaN(eeg.nbchan, eeg.nbchan);
            pMatrix = NaN(eeg.nbchan, eeg.nbchan);
            for chan1 = 1:eeg.nbchan
                data1 = eeg.data(chan1, startIndex1:endIndex1)';
                %     for chan2 = (chan1+1):eeg.nbchan
                for chan2 = 1:eeg.nbchan
                    startIndex2 = searchRange + startIndex1;
                    endIndex2 = searchRange + endIndex1;
                    remove = startIndex2 < 1;
                    remove = remove | (endIndex2 > eeg.pnts);
                    if(any(remove))
                        startIndex2(remove) = [];
                        endIndex2(remove) = [];
                    end
                    bestP = realmax;
                    bestRho = NaN;
                    bestLagFrames = NaN;
                    for searchIndex = 1:length(startIndex2)
                        startI2 = startIndex2(searchIndex);
                        endI2 = endIndex2(searchIndex);
                        data2 = eeg.data(chan2, startI2:endI2)';
                        [rho, p] = corr(data1, data2);
                        if(p < bestP)
                            lagFrames = startI2 - startIndex1;
                            bestP = p;
                            bestRho = rho;
                            bestLagFrames = lagFrames;
                        end
                    end
                    lagMatrix(chan1, chan2) = bestLagFrames;
                    pMatrix(chan1, chan2) = bestP;
                    rhoMatrix(chan1, chan2) = bestRho;
                end
            end
            summary.lagMatrix = lagMatrix;
            summary.pMatrix = pMatrix;
            summary.rhoMatrix = rhoMatrix;
            summary.parameters.inputFile = filename;
            summary.parameters.sampleRate = eeg.srate;
            summary.parameters.windowStartFrame = startIndex1;
            summary.parameters.windowEndFrame = endIndex1;
            summary.parameters.channelLocations = eeg.chanlocs;
            summary.parameters.searchLimits = searchLimits;
            summary.parameters.searchFrameInterval = searchFrameInterval;
            
            save(outputFile, 'summary', '-v7.3')
        end
    end
    %     end
end

figure;
imagesc(lagMatrix);
title('lag');
figure;
imagesc(rhoMatrix);
title('rho');
figure;
imagesc(pMatrix);
title('p');

[clusterRho, perm] = clusterMatrix(rhoMatrix);
permLab = {eeg.chanlocs.labels};
permLab = permLab(perm);

figure;
imagesc(clusterRho);
title('rho (clustered)');
set(gca, 'xtick', 1:eeg.nbchan);
set(gca, 'xticklabel', permLab);
set(gca, 'ytick', 1:eeg.nbchan);
set(gca, 'yticklabel', permLab);


clusterLag = clusterMatrix(lagMatrix, perm);
figure;
imagesc(clusterLag);
title('lag (clustered)');
set(gca, 'xtick', 1:eeg.nbchan);
set(gca, 'xticklabel', permLab);
set(gca, 'ytick', 1:eeg.nbchan);
set(gca, 'yticklabel', permLab);


% tab=tabulate(diff(positions));
% tab(tab(:,2)==0,:) = [];
% tab

