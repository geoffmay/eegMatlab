function trainingData = getConvolvedEegNeuralNetworkData(edfFilename)
%compute phase slope info
% edfFilename = 'GeoffTestEEG2-edf.edf';
%1: DMN
%2: Visual
%3: FP
%5: DAN
%7: VAN/Language
%8: Salience
%9: Cingulo-Opercular
%10: Motor hand
%11: Motor mouth
%12: Aud
%15: PMN
%16: CAN

doPhase = false;

% edfFilename = 'GeoffTestEEG2-edf.edf'

keepNetworks = [1 3];

restNetworkTimecourse = textread('C:\Vision\Raw Files\Geoff EEG test\mri\sub-GeoffEEGTest_RSFC_network_timecourses.txt');
restNetworkIndices = restNetworkTimecourse(1,:);
restNetworkTimecourse = restNetworkTimecourse(2:end,keepNetworks);


metadata = readCsv('C:\Vision\Raw Files\Geoff EEG test\mri\sub-GeoffEEGTest_metadata.csv');



markerFilename = strrep(edfFilename, '-edf.edf', '_Pulse Artifact Correction.Markers');
locations = fileLocations;
markers = loadVolumeMarkers(fullfile(locations.brainVision, markerFilename));
if(doPhase)
    phasePath = ['C:\Vision\Raw Files\Geoff EEG test\history\', edfFilename, 'phase.mat'];
    phase = load(phasePath);
    convPhase = phase.convPhase(markers(1:size(restNetworkTimecourse, 1)), :);
end

[eeg, dropFrames] = loadBrainvisionEdf(fullfile(locations.brainVision, edfFilename));

%load convolved eeg signals.
convolvedFilename = ['C:\Users\Neuro\Documents\MATLAB\processed\boldEeg\convolved\', edfFilename, 'convolved.mat'];
if(~exist(convolvedFilename, 'file'))
    % tic;convolvedEeg = convolveEegAtVolumeMarkers(edfFilename, eeg);toc;
    % tic;convolvedEeg = convolveEegMeasures(edfFilename);toc;
    %tic;convolvedEeg = consolidateDerivedParameters(edfFilename);toc;
    tic;convolvedEeg = loadSplitConvolution(edfFilename);toc;
    save(convolvedFilename, 'convolvedEeg', '-v7.3');
else
    tic;load(convolvedFilename);toc;
end

%drop volumes where MRI gradient artifact could affect convolved signal
keepVolumes = false(size(convolvedEeg.timeSeconds));
keepVolumes(1:size(restNetworkTimecourse,1)) = true;
rises = find(diff(dropFrames) == 1);
falls = find(diff(dropFrames) == -1);
if(falls(1) < rises(1))
    rises = [1, rises];
end
if(length(rises) > length(falls))
    falls = [falls, eeg.pnts];
end
rises = rises ./ eeg.srate;
falls = falls ./ eeg.srate;
for i = 1:length(rises)
    dropped = convolvedEeg.timeSeconds > rises(i) & convolvedEeg.timeSeconds < falls(i);
    keepVolumes(dropped) = false;
end

%target alpha power/coherence only
freqTarget = '9-12';
% measureTarget = 'ower';
measureTarget = '';

if(length(freqTarget) > 0)
    isFreq = cellfun(@length, strfind(convolvedEeg.labels, freqTarget)) > 0;
else
    isFreq = true(size(convolvedEeg.labels));
end
if(length(measureTarget) > 0)
    isMeas = cellfun(@length, strfind(convolvedEeg.labels, measureTarget)) > 0;
else
    isMeas = true(size(convolvedEeg.labels));
end
desiredEeg = isFreq & isMeas;
eegInput = convolvedEeg.matrix(keepVolumes, desiredEeg);
labels = convolvedEeg.labels(desiredEeg);

if(doPhase)
    eegInput = [eegInput, convPhase];
end

%drop bad MRI volumes (eeg was already dropped)
restNetworkTimecourse = restNetworkTimecourse(keepVolumes(1:size(restNetworkTimecourse,1)) ,:);

trainingData.eegInput = eegInput;
trainingData.mriOutput = restNetworkTimecourse;
trainingData.eegColumnLabels = labels;


if(false)
    %pca
    [eegPca.COEFF, eegPca.SCORE, eegPca.LATENT, eegPca.TSQUARED, eegPca.EXPLAINED, eegPca.MU] = pca(eegInput);
    eegPcaInput = eegPca.SCORE(:, 1:100);
    
    %ica
    [eegIca.weights,eegIca.sphere,eegIca.compvars,eegIca.bias,eegIca.signs,eegIca.lrates,eegIca.activations] = runica(eegInput');
end