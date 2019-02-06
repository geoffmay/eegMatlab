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


testOnly = true;
doPhase = false;

edfFilename = 'GeoffTestEEG2-edf.edf'

keepNetworks = [1 3];

restNetworkTimecourse = textread('C:\Vision\Raw Files\Geoff EEG test\mri\sub-GeoffEEGTest_RSFC_network_timecourses.txt');
restNetworkIndices = restNetworkTimecourse(1,:);
restNetworkTimecourse = restNetworkTimecourse(2:end,keepNetworks);
%restNetworkTimecourse(:, all(isnan(restNetworkTimecourse))) = [];


metadata = readCsv('C:\Vision\Raw Files\Geoff EEG test\mri\sub-GeoffEEGTest_metadata.csv');



markerFilename = strrep(edfFilename, '-edf.edf', '_Pulse Artifact Correction.Markers');
locations = fileLocations;
markers = loadVolumeMarkers(fullfile(locations.brainVision, markerFilename));
if(doPhase)
    phasePath = ['C:\Vision\Raw Files\Geoff EEG test\history\', edfFilename, 'phase.mat'];
    phase = load(phasePath);
    convPhase = phase.convPhase(markers(1:size(restNetworkTimecourse, 1)), :);
end
eeg = loadBrainvisionEdf(fullfile(locations.brainVision, edfFilename));

if(testOnly)
    pointCount = 10000;
    eeg.pnts = pointCount;
    eeg.data(:, (pointCount+1):end) = [];
    eeg.times(pointCount+1) = [];
end

%load convolved eeg signals.
convolvedFilename = ['C:\Users\Neuro\Documents\MATLAB\processed\MrEegLink\', edfFilename, 'convolved.mat'];
if(~exist(convolvedFilename, 'file'))
    
    tic;eegToConvolution(edfFilename);toc;
    if(testOnly)
        convolvedEeg
    end

    %tic;convolvedEeg = consolidateDerivedParameters(edfFilename);toc;
    tic;convolvedEeg = loadSplitConvolution(edfFilename);toc;
    save(convolvedFilename, 'convolvedEeg', '-v7.3');
else
    tic;load(convolvedFilename);toc;
end

target = '9-12';
rowKeep = 1:size(restNetworkTimecourse,1);
if(length(target) > 0)
    isHit = cellfun(@length, strfind(convolvedEeg.labels, target)) > 0;
    eegInput = convolvedEeg.matrix(rowKeep, isHit);
else
    eegInput = convolvedEeg.matrix(rowKeep, :);
end
if(doPhase)
    eegInput = [eegInput, convPhase];
end



if(false)
    %pca
    [eegPca.COEFF, eegPca.SCORE, eegPca.LATENT, eegPca.TSQUARED, eegPca.EXPLAINED, eegPca.MU] = pca(eegInput);
    eegPcaInput = eegPca.SCORE(:, 1:100);
    
    %ica
    [eegIca.weights,eegIca.sphere,eegIca.compvars,eegIca.bias,eegIca.signs,eegIca.lrates,eegIca.activations] = runica(eegInput');
end