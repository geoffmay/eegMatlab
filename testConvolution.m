

if(~exist('corrLagPca', 'var'))
corrLagFolder = 'C:\Vision\Raw Files\Geoff EEG test\export\corrLag';
corrLagFiles = dir(corrLagFolder);
corrLagFiles([corrLagFiles.isdir]) = [];
corrLagCounter = 0;
corrLagTimeFlatCourse = [];
corrLagStartTimes = NaN(1, length(corrLagFiles));
corrLagEndTimes = NaN(1, length(corrLagFiles));
for corrLagFileCounter = 1:length(corrLagFiles)
    fprintf('\ncorrLag file %d of %d', corrLagFileCounter, length(corrLagFiles));
    corrLagData = load(fullfile(corrLagFiles(corrLagFileCounter).folder, corrLagFiles(corrLagFileCounter).name));
    if(isfield(corrLagData, 'summary'))
        corrLagFlat = reshape(corrLagData.summary.lagMatrix, [1, numel(corrLagData.summary.lagMatrix)]);
        if(corrLagFileCounter == 1)
            corrLagTimeFlatCourse = NaN(length(corrLagFiles), size(corrLagFlat, 2));
        end
        corrLagStartTimes(corrLagFileCounter) = corrLagData.summary.parameters.windowStartFrame;
        corrLagTimeFlatCourse(corrLagFileCounter, :) = corrLagFlat;
        corrLagStartTimes(corrLagFileCounter) = corrLagData.summary.parameters.windowStartFrame / corrLagData.summary.parameters.sampleRate;
        corrLagEndTimes(corrLagFileCounter) = corrLagData.summary.parameters.windowEndFrame / corrLagData.summary.parameters.sampleRate;
        chanlocs = corrLagData.summary.parameters.channelLocations;
    end
end
remove = isnan(corrLagTimeFlatCourse(:,1));
corrLagTimeFlatCourse(remove, :) = [];
corrLagStartTimes(remove) = [];
corrLagEndTimes(remove) = [];


[corrLagPca.COEFF, corrLagPca.SCORE, corrLagPca.LATENT, corrLagPca.TSQUARED, corrLagPca.EXPLAINED, corrLagPca.MU] = pca(corrLagTimeFlatCourse);
end

for i = 1:20;
    compFlat = corrLagPca.COEFF(:, i);
    comp = reshape(compFlat, [sqrt(length(compFlat)), sqrt(length(compFlat))]);
    
    if(false)
        close all;
        [clustComp, clustPerm] = clusterMatrix(comp);
        clustChan = corrLagData.summary.parameters.channelLocations(clustPerm);
        figure;
        imagesc(clustComp);
        set(gca, 'xtick', 1:length(clustPerm));
        set(gca, 'xticklabel', {clustChan.labels});
        set(gca, 'ytick', 1:length(clustPerm));
        set(gca, 'yticklabel', {clustChan.labels});
        mytitle = sprintf('max correlatopm time lag component %d (%f%% explained)', i, corrLagPca.EXPLAINED(i));
        
        figure;
        % plot(corrLagPca.SCORE(:, i));
        compBold = convolveHrf(corrLagPca.SCORE(:, i), 4);
        compBold(length(corrLagStartTimes)+1:end) = [];
        hold on;
        plot(corrLagStartTimes, compBold);
        mytitle = sprintf('max correlatopm time lag component %d (%f%% explained)', i, corrLagPca.EXPLAINED(i));
        title(mytitle);
    end
    
    figure;
    topoplot(sum(comp, 1), chanlocs);
    mytitle = sprintf('max correlatopm time lag component %d (%f%% explained)', i, corrLagPca.EXPLAINED(i));
    title(mytitle);
end
tilefigs;

tic;
convComp2 = convolveHrf(corrLagPca.SCORE(:,2), 4);
toc;

convEeg = convComp2;
volumeMarkerFrames = loadVolumeMarkers('C:\Vision\Raw Files\Geoff EEG test\export\GeoffTestEEG2_Pulse Artifact Correction.Markers');
volumeMarkerTimes = volumeMarkerFrames / 250;
% volumeMarkersDownsampled = volumeMarkers .* (250 / 5000);
convolvedEegAtVolumeMarkers = interp1(corrLagStartTimes, convEeg(1:length(corrLagStartTimes)), volumeMarkerTimes)';
networkAvgs = textread('C:\Vision\Raw Files\Geoff EEG test\export\networkavgs.txt');
keep = all(~isnan(networkAvgs), 2);

for i = 1:size(networkAvgs, 2)
    [rhos(i), ps(i)] = corr(convolvedEegAtVolumeMarkers(keep), networkAvgs(keep, i));
end


%do basic stuff

eeg = loadBrainvisionEdf('C:\Vision\Raw Files\Geoff EEG test\export\GeoffTestEEG2-edf.edf');
surfCoh = deriveRobiCoherenceMatrix(eeg);
saveMemory = true;
if(saveMemory)
    clear corrLagPca eeg corrLagTimeFlatCourse
    dsRatio = 5;
    for i = 1:length(surfCoh.coh)
        fprintf('\ncoh %d of %d', i, length(surfCoh.coh));
        temp = downsample(surfCoh.coh(i).coherence, dsRatio, 1);
        surfCoh.coh(i).coherence = temp;
    end
    for i = 1:length(surfCoh.channels)
        fprintf('\nchannel %d of %d', i, length(surfCoh.channels));
        surfCoh.channels(i).absolutePower = downsample(surfCoh.channels(i).absolutePower, dsRatio, 1);
        surfCoh.channels(i).relativePower = downsample(surfCoh.channels(i).relativePower, dsRatio, 1);
    end
    surfCoh.x = downsample(surfCoh.x, dsRatio, 1);
end

save('C:\Vision\Raw Files\eegtest\workbench2.mat', '-v7.3');
cohMat = convertCoherenceStructToMatrix(surfCoh.coh, surfCoh.freqInfo, surfCoh.channels);



save('C:\Vision\Raw Files\eegtest\workbench2.mat', '-v7.3');


%don't do phase stuff yet
if(false)

if(false)
    x1 = [2 6 4 9];
    n1 = [2 3 4 5];
    x2 = [2 7 5];
    n2 = [8 9 10];
    
    [y,n] = convolution(x1,n1,x2,n2);% will give the following output-
    % y = [4 26 60 76 83 45];
    % n = [10 11 12 13 14 15];
    
    [yb,nb] = convolution(x2,n2,x1,n1);% will give the following output-
end

if(false)
    x = zeros(1, 10000);
    x(50:100) = 1;
    x(8000:9000) = 1;
    
    [xBold, xTimes] = convolveHrf(x, 250);
    close all;
    plot(xTimes, xBold);
end

phaseFolder = 'C:\Vision\Raw Files\eegtest\export\phaseSlope';
phaseFiles = dir(phaseFolder);
phaseFiles([phaseFiles.isdir]) = [];
timeLags = [];
timeLagCounter = 0;
for phaseFileCounter = 1:length(phaseFiles)
    fprintf('\nphase file %d of %d', phaseFileCounter, length(phaseFiles));
    phaseData = load(fullfile(phaseFiles(phaseFileCounter).folder, phaseFiles(phaseFileCounter).name));
    if(isfield(phaseData, 'topography'))
        timeLag = phaseData.topography.estimatedTimeLag;
        remove = isnan(timeLag(:,1,1));
        timeLag(remove, :, :) = [];
        if(phaseFileCounter == 1)
            timeLags = NaN(size(timeLag,1) * length(phaseFiles), size(timeLag, 2), size(timeLag,3));
        end
        timeLags((timeLagCounter+1):(timeLagCounter+size(timeLag,1)), :, :) = timeLag;
        timeLagCounter = timeLagCounter + size(timeLag,1);
    end
end
timeLags(timeLagCounter+1:end, :, :) = [];
[phasePca.COEFF, phasePca.SCORE, phasePca.LATENT, phasePca.TSQUARED, phasePca.EXPLAINED, phasePca.MU] = pca(timeLags(:, :, 1));

i = 2;
component = phasePca.COEFF(:, i);
figure;
topoplot(component, topography.chanlocs);
title(sprintf('phase lag pca component %d (%f%% explained)', i, phasePca.EXPLAINED(i)));
brainActivity = phasePca.SCORE(:, i);
brainBold = convolveHrf(brainActivity, 250);
figure;
plot(brainActivity);
end

