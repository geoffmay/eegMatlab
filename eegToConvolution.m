function convCoh = eegToConvolution(edfFilename)

locations = fileLocations;
if(~isstruct(edfFilename))
    ind = strfind(edfFilename, 'EEG_data_sub-');
    if(length(ind) > 0)
        numText = edfFilename((ind(end) + length('EEG_data_sub-')):end);
        underscore = strfind(numText, '_');
        numText(underscore:end) = [];
        folder = fullfile(locations.ghermanFolder, ['sub-GhermanPhiliastides', numText], 'EEG');
        eeg = loadGhermanEeg(fullfile(folder, edfFilename));
    else
        locations = fileLocations();
        folder = locations.brainVision;
        eeg = loadBrainvisionEdf(fullfile(folder, edfFilename));
    end
else
    eeg = edfFilename;
    edfFilename = eeg.filepath;
end


networkFilename = strrep(eeg.filename, '-edf.edf', '_network_timecourses.txt');
networkActivity = textread(['C:\Vision\Raw Files\Geoff EEG test\mri\' networkFilename]);
networkLabels = networkActivity(1, :);
networkActivity = networkActivity(2:end, :);

targetNetworks = networkActivity(:, [1 3]);


%power and coherence measures
[cohInfo, summary] = deriveRobiCoherenceMatrix(eeg);
saveSplitCohFile(cohInfo);
clear cohInfo summary;
convCoh = convolveAtVolumeMarkers(edfFilename);



if(false)
    %bandpass filter
    %this makes the eeg data resemble gherman et al
    eeg.data = oldData;
    for i = 1:size(eeg.data,1)
        fprintf('channel %d\n', i);
        sig = eeg.data(i,:);
        %     filtSig = bandpass(sig, [1, 40], 1000);
        filtSig = bandpass(sig, [1, 40], 250);
        eeg.data(i,:) = filtSig;
    end
    
    %ica
    %experimental
    [eegIca.weights,eegIca.sphere,eegIca.compvars,eegIca.bias,eegIca.signs,eegIca.lrates,eegIca.activations] = runica(eeg.data);
    icaEeg = eeg;
    icaEeg.data = eegIca.activations;
    
    %convolve
    convovled = convolveEegAtVolumeMarkers(eeg.filename, icaEeg);
    icaInput = convovled.signal(:, 1:(132*3))';

end

